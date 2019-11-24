###################################################
## scenarios.jl
##      Compute routing solutions for a single school
## Authors: Arthur Delarue, Sébastien Martin, 2018
###################################################

"""
    The main euclidian traveltime function.
    Select the time that is the closest to `time` and return the euclidian travel-time
    for the corresponding speed.
    o:: Origin
    d:: Destination
"""
 # using JuMP, Gurobi
using JuMP, Cbc
include("problem.jl")
include("load.jl")
#Sinclude("SchoolBusRouting.jl")
import Statistics: mean
using Base:GC
using ProgressMeter
using Random


data = loadSyntheticBenchmark("data/input/CSCB01/Schools.txt",
                                  "data/input/CSCB01/Stops.txt")

# data = loadSyntheticBenchmark("../data/input/CSCB01/Schools.txt",
                                  # "../data/input/CSCB01/Stops.txt")
data.stops[1]

function traveltime(data::SchoolBusData, o::Point, d::Point)
    if data.params.metric == MANHATTAN
        return manhattandistance(o,d) / data.params.velocity
    else
        return euclideandistance(o,d) / data.params.velocity
    end
end



"""
    Manhattan Distance between two points, in meters
"""
function manhattandistance(o::Point, d::Point)
    return abs(o.x - d.x) + abs(o.y - d.y)
end

"""
    Euclidean distance between two points, in meters
"""
function euclideandistance(o::Point, d::Point)
    return sqrt((o.x - d.x) ^ 2 + (o.y - d.y) ^ 2)
end

"""
School, Yard, Stop  2^3
"""
traveltime(data::SchoolBusData, o::Yard,   d:: Stop)   = traveltime(data, o.position, d.position)
traveltime(data::SchoolBusData, o::School, d:: Stop)   = traveltime(data, o.position, d.position)
traveltime(data::SchoolBusData, o::Stop,   d:: Stop)   = traveltime(data, o.position, d.position)
traveltime(data::SchoolBusData, o::Stop,   d:: School) = traveltime(data, o.position, d.position)
traveltime(data::SchoolBusData, o::School, d:: Yard)   = traveltime(data, o.position, d.position)
traveltime(data::SchoolBusData, o::Yard,   d:: School) = traveltime(data, o.position, d.position)
traveltime(data::SchoolBusData, o::Stop,   d:: Yard)   = traveltime(data, o.position, d.position)

"""
    Get number of students at the stop
"""
nStudents(data::SchoolBusData, stop::Stop) = stop.nStudents
nStudents(data::SchoolBusData, school::Int, stop::Int) = data.stops[school][stop].nStudents
#
# bus_data = loadSchoolsReduced("../data/input/CSCB01/Schools.txt")
#
# bus_data[1].position
data.schools[1].position
"""
    Get time that bus must remain at stop
"""
function stopTime(data::SchoolBusData, stop::Stop)
    return data.params.constant_stop_time + data.params.stop_time_per_student *
    stop.nStudents
end

# stop_data = CSV.read("../data/input/CSCB01/Stops.txt")
# stops = stop_data[1,:]
# stops = data.stops[1][1]
# stopTime(data, stops)
"""
    Get the maximum allowed travel time for a given bus stop
"""
function maxTravelTime(data::SchoolBusData, stop::Stop)
    return data.params.max_time_on_bus
end

maxTravelTime(data, data.stops[2][3]) # constant time 2700s
"""
    The current state of the greedy algorithm
    Contains:
        - the current route
        - an occupancy object that tells us how many students are currently on the bus
        - for each stop, the amount of extra time we can leave the students on the bus
        - for each stop, the current time spent on the bus
        - the total service time of the route (including time of first stop)
"""
struct GreedyState
    "the current route"
    route::Route
    "the current occupancy"
    nStudents::Int
    "amount of extra time students at each stop can spend on bus"
    slackTimes::Vector{Float64}
    "amount of time students at each stop are already spending on bus"
    stopTimes::Vector{Float64}
    "total time of current route"
    routeTime::Float64
end


"""
    Greedy single-school routing heuristic
    Args:
        - the data
        - the ID of the school
        - the maximum time allowed on a route
        - the penalty (in seconds) for adding D2D students to routes
    Returns:
        - all routes for the school as a Vector{Route}
"""

maxRouteTime = 2700.0

Random.seed!(1234)

function greedy(data::SchoolBusData, schoolID::Int, maxRouteTime::Float64)
    routes = Route[]
    availableStops = trues(length(data.stops[schoolID]))
    while sum(availableStops) > 0
        # randomly select a starting stop
        startingStopID = rand((1:length(availableStops))[availableStops])
        # create initial route (select initial bus)
        currentState = initialRoute(data, schoolID, startingStopID,
                                    length(routes)+1)
        availableStops[startingStopID] = false
        while true
            bestStopID = 0
            bestInsertId = -1
            bestTimeDiff = Inf
            for stopID in collect(1:length(availableStops))[availableStops]
                # check the students number when adding a stop
                if (data.stops[schoolID][stopID].nStudents + currentState.nStudents <=
                    data.params.bus_capacity)

                    insertId, timeDiff = bestInsertion(data, schoolID, stopID,
                                                       currentState, maxRouteTime)
                    if timeDiff < bestTimeDiff
                        bestTimeDiff = timeDiff
                        bestStopID = stopID
                        bestInsertId = insertId
                    end
                end
            end

            if bestTimeDiff < Inf
                currentState = buildRoute(data, currentState, schoolID,
                                          bestStopID, bestInsertId)
                availableStops[bestStopID] = false
            else
                push!(routes, currentState.route)
                break
            end
        end
    end
    return routes
end

schoolID = 1
stopID = 1

# maxRouteTime = 2700.0
# all_routes = greedy(data, 1, maxRouteTime)
# all_routes[1]



"""
    Get total travel time difference when adding stop to the route, together with
        a "cost" of how good this insertion is.
    The cost is equal to the total amount of time by which this insertion
        increases the length of the route, not including stop time
    Args:
        - data          : the data object, assumes it has the stops computed
        - schoolID      : the school we are routing
        - newStopID     : the stop we are trying to insert
        - cs            : the current GreedyState
        - maxRouteTime  : the maximal allowed length of the route
    Returns:
        - insertID      : where we would optimally insert it
        - bestTimeDiff  : the optimal time difference imposed on the route
   Comment:
       when adding a new stop, if we should add in the 1st,
       or in the middle or the last stop before school.
"""
function bestInsertion(data::SchoolBusData,
                       schoolID::Int,
                       newStopID::Int,
                       cs::GreedyState,
                       maxRouteTime::Float64)
    newStop      = data.stops[schoolID][newStopID]
    school       = data.schools[schoolID]
    bestTimeDiff = Inf
    insertId     = -1
    # before first stop
    timeDiff = traveltime(data, newStop, data.stops[schoolID][cs.route.stops[1]])
    totTime = timeDiff + cs.routeTime + stopTime(data, newStop)
    if (totTime - stopTime(data, newStop) <= maxTravelTime(data, newStop) &&
                totTime <= maxRouteTime)
        bestTimeDiff = timeDiff
        insertId = 0
    end
    # between stops
    # loop over all the stops in the middle
    # local search to get the best feasible route.

    for i = 1:(length(cs.route.stops)-1)
        timeToNextStop = traveltime(data, newStop, data.stops[schoolID][cs.route.stops[i+1]])
        timeDiff = traveltime(data, data.stops[schoolID][cs.route.stops[i]], newStop) +
                   timeToNextStop -
                   traveltime(data, data.stops[schoolID][cs.route.stops[i]],
                      data.stops[schoolID][cs.route.stops[i+1]])
        if timeDiff < bestTimeDiff
            totTime = timeDiff + cs.routeTime + stopTime(data, newStop)
            # check feasibility for all stops preceding this one
            isFeasible = (timeDiff + stopTime(data, newStop) <= cs.slackTimes[i])
            # check feasibility for the potential new stop
            isFeasible = isFeasible && (cs.stopTimes[i+1] + timeToNextStop +
                                        stopTime(data,
                                             data.stops[schoolID][cs.route.stops[i+1]]) <=
                                        maxTravelTime(data, newStop))
            isFeasible = isFeasible && totTime <= maxRouteTime
            if isFeasible
                bestTimeDiff = timeDiff
                insertId = i
            end
        end
    end
    # after last stop
    timeDiff =  traveltime(data, data.stops[schoolID][cs.route.stops[end]], newStop) +
                traveltime(data, newStop, school) -
                traveltime(data, data.stops[schoolID][cs.route.stops[end]], school)
    if timeDiff < bestTimeDiff
        totTime = timeDiff + cs.routeTime + stopTime(data, newStop)
        isFeasible = timeDiff + stopTime(data, newStop) <= cs.slackTimes[end] &&
                        (traveltime(data, newStop, school) <= maxTravelTime(data, newStop))
        isFeasible = isFeasible && totTime <= maxRouteTime
        if isFeasible
            bestTimeDiff = timeDiff
            insertId = length(cs.route.stops)
        end
    end

    return insertId, bestTimeDiff
end

"""
    Given new stop ID and current route/GreedyState, update route
    Args:
        - data          : the data object, assumes it has the stops computed
        - cs            : the current GreedyState
        - schoolID      : the school we are routing
        - newStopID     : the stop we are inserting
        - insertID      : where we insert it
    Returns:
        - a new GreedyState with the updated route
"""
function buildRoute(data::SchoolBusData,
                    cs::GreedyState,
                    schoolID::Int,
                    newStopID::Int,
                    insertID::Int)
    newStop = data.stops[schoolID][newStopID]
    school = data.schools[schoolID]
    if insertID == 0 # at the beginning
        nextStop = data.stops[schoolID][cs.route.stops[1]]
        "amount of time students at each stop are already spending on bus"
        newStopTimeOnBus = traveltime(data, newStop, nextStop) +
                           cs.stopTimes[1] + stopTime(data,nextStop)
        timeDiff = 0.
    elseif insertID == length(cs.route.stops) # at the end
        previousStop = data.stops[schoolID][cs.route.stops[end]] # last stop on route
        newStopTimeOnBus = traveltime(data, newStop, school) # time from new stop to school
        timeDiff = newStopTimeOnBus +
                   traveltime(data, previousStop, newStop) -
                   traveltime(data, previousStop, school)
    else # in the middle
        previousStop = data.stops[schoolID][cs.route.stops[insertID]]
        nextStop = data.stops[schoolID][cs.route.stops[insertID+1]]
        timeToNextStop = traveltime(data, newStop, nextStop)
        timeDiff = traveltime(data, previousStop, newStop) +
                   timeToNextStop -
                   traveltime(data, previousStop, nextStop)
        newStopTimeOnBus = timeToNextStop + cs.stopTimes[insertID+1] + stopTime(data, nextStop)
    end

    # it inclueds all stops' stoptime and slackTimes
    # if insertId == 0, then 1:insertID will aalwyas return 0.
    newStopTimes = vcat(cs.stopTimes[1:insertID] .+ stopTime(data, newStop) .+ timeDiff,
                        [newStopTimeOnBus],
                        cs.stopTimes[(insertID+1):end])
    newStops = vcat(cs.route.stops[1:insertID], [newStopID], cs.route.stops[(insertID+1):end])
    newSlackTimes = vcat(cs.slackTimes[1:insertID] .- stopTime(data, newStop) .- timeDiff,
                         [maxTravelTime(data, newStop)] .- newStopTimeOnBus,
                         cs.slackTimes[(insertID+1):end])
    # fix slack time property
    for i=eachindex(newSlackTimes)
        if i > 1
            newSlackTimes[i] = min(newSlackTimes[i-1], newSlackTimes[i])
        end
    end
    routeTime = newStopTimes[1] + stopTime(data, data.stops[schoolID][newStops[1]])
    return GreedyState(Route(cs.route.id, newStops),
                       cs.nStudents + newStop.nStudents,
                       newSlackTimes, newStopTimes, routeTime)
end

# cs = initial_route_test
# println(cs)

# newStopID = 2
# insertID = 0
# schoolID = 1
# route_example = buildRoute(data, cs,
#  schoolID, newStopID, insertID )
# println(route_example)
# cs = route_example
#
# newStopID = 2
# insertID = 1
# schoolID = 1
#
# route_example1 = buildRoute(data, route_example,
#  schoolID, newStopID, insertID )
#  println(route_example1)

"""
    Create initial route for greedy
    Args:
        - data          : the data object, assumes it has the stops computed
        - schoolID      : the school we are routing
        - stopID        : the stop we start with
        - routeID       : the ID of the route we are creating
    Returns a GreedyState object
"""
schoolID = 1
stopID = 1
function initialRoute(data::SchoolBusData,
                      schoolID::Int,
                      stopID::Int,
                      routeID::Int)
    timeOnBus = traveltime(data, data.stops[schoolID][stopID], data.schools[schoolID])
    nStudents = data.stops[schoolID][stopID].nStudents
    slackTimes = [maxTravelTime(data, data.stops[schoolID][stopID]) - timeOnBus]
    "amount of time students at each stop are already spending on bus"
    stopTimes = [timeOnBus]
        routeTime = timeOnBus + stopTime(data, data.stops[schoolID][stopID])
    return GreedyState(Route(routeID, [stopID]), nStudents, slackTimes, stopTimes, routeTime)
end


initial_route_test = initialRoute(data,data.stops[1][1].id, data.stops[1][1].id, 1)
println(initial_route_test)
"""
    Returns the sum of each student's travel times on this route
"""

function sumIndividualTravelTimes(data::SchoolBusData, schoolID::Int, r::Route)
    allStops = data.stops[schoolID]
    numStops = length(r.stops)
    t = traveltime(data, allStops[r.stops[numStops]], data.schools[schoolID]) * numStops
    while numStops > 1
        t += (traveltime(data, allStops[r.stops[numStops-1]], allStops[r.stops[numStops]]) +
              stopTime(data, allStops[r.stops[numStops]])) * (numStops - 1)
        numStops -= 1
    end
    return t
end
#
# r = all_routes[2]
# sumIndividualTravelTimes(data, 1, r)
function sumIndividualTravelTimes(data::SchoolBusData, schoolID::Int,
                                  routes::Vector{Route})
    return sum(sumIndividualTravelTimes(data, schoolID, route) for route in routes)
end

# sumIndividualTravelTimes(data, 1, all_routes)


"""
    Returns the total travel time from the time the first student enters the bus to the time of pickup/dropoff in school
"""
function serviceTime(data::SchoolBusData, schoolID::Int, stoplist::Vector{Int})
    allStops = data.stops[schoolID]
    numStops = length(stoplist)
    t = traveltime(data, allStops[stoplist[numStops]], data.schools[schoolID])
    while numStops > 1
        t += (traveltime(data, allStops[stoplist[numStops-1]], allStops[stoplist[numStops]]) +
              stopTime(data, allStops[stoplist[numStops]]))
        numStops -= 1
    end
    t += stopTime(data, allStops[stoplist[1]])
    return t
end
serviceTime(data::SchoolBusData, schoolID::Int, r::Route) = serviceTime(data, schoolID, r.stops)

# serviceTime(data, 1, r)



"""
    Simple route representation : just list of IDs and associated cost
"""
struct FeasibleRoute
    "The list of stop ids"
    stopIds::Vector{Int}
    "The cost associated with the route"
    cost::Float64
end

function FeasibleRoute(data::SchoolBusData, schoolID::Int, r::Route)
    return FeasibleRoute(r.stops, sumIndividualTravelTimes(data, schoolID, r))
end

# a_feasible_route = FeasibleRoute(data, 1, r)

"""
    Stores a list of routes in a way that makes column generation easy
"""
struct FeasibleRouteSet
    "The list of FeasibleRoutes available"
    list::Vector{FeasibleRoute}
    "The set of stopIds list (one for each route)"
    set::Set{Vector{Int}}
    "For each stop, the index of the routes in the list that go through this stop"
    atStop::Vector{Vector{Int}}

    FeasibleRouteSet(data::SchoolBusData, schoolID::Int) =
        new(FeasibleRoute[], Set{Vector{Int}}(), [Int[] for i = 1:length(data.stops[schoolID])])
end

# a_feasible_set = FeasibleRouteSet(data, 1)
"""
    Generate N random greedy solutions, and combines them smartly to get the best
    possible solution.
"""
function greedyCombined(data::SchoolBusData,
                        schoolID::Int,
                        N::Int,
                        maxRouteTimeLower::Float64,
                        maxRouteTimeUpper::Float64,
                        λ::Float64;
                        # ,
                        # env::Gurobi.Env=Gurobi.Env();
                        args...)
    routeList = generateRoutes(data, schoolID, N, maxRouteTimeLower, maxRouteTimeUpper)
    routes = FeasibleRouteSet(data, schoolID)
    addRoute!(routes, routeList)
    # selectedRoutes = bestRoutes(data, schoolID, routes, λ, env; args...)
    selectedRoutes = bestRoutes(data, schoolID, routes, λ;args... )

    return buildSolution(data, schoolID, routes, selectedRoutes)
end

# greedyCombined(data, schoolID, N, maxRouteTimeLower, maxRouteTimeUpper, λ  )

function greedyCombined(data::SchoolBusData,
                        schoolID::Int,
                        startRoutes::Vector{Route},
                        N::Int,
                        maxRouteTimeLower::Float64,
                        maxRouteTimeUpper::Float64,
                        λ::Float64;
                        # change it back if gurobi is available
                        # ,
                        # env::Gurobi.Env=Gurobi.Env();
                        args...
                        )
    routeList = generateRoutes(data, schoolID, N, maxRouteTimeLower, maxRouteTimeUpper)
    routes = FeasibleRouteSet(data, schoolID)
    addRoute!(routes, routeList)
    addRoute!(routes, collect(FeasibleRoute(data, schoolID, r) for r in startRoutes))
    # selectedRoutes = bestRoutes(data, schoolID, routes, λ, env; args...)
    selectedRoutes = bestRoutes(data, schoolID, routes, λ;args... )

    return buildSolution(data, schoolID, routes, selectedRoutes)
end


# startRoutes =  greedy(data, schoolID, Inf)

# greedyCombined(data, schoolID, startRoutes, N, maxRouteTimeLower, maxRouteTimeUpper, λ  )

"""
    Same as greedy combined, but iterates it to keep improving the best solution
"""
function greedyCombinedIterated(data::SchoolBusData,
                                schoolID::Int,
                                maxRouteTimeLower::Float64,
                                maxRouteTimeUpper::Float64,
                                nGreedy::Int,
                                nIteration::Int,
                                λ::Float64;
                                verbose::Bool=false,
                                # ,
                                args...
                                )
    # env = Gurobi.Env()
    routes = greedy(data, schoolID, Inf)
    verbose && @printf("Iteration 0: %d buses, %2.fs\n", length(routes),
                        sumIndividualTravelTimes(data, schoolID, routes))
    for i=1:nIteration
        routes = greedyCombined(data, schoolID, routes, nGreedy,
                                maxRouteTimeLower, maxRouteTimeUpper, λ;
                                # , env;
                                args...
                                )

        verbose && @printf("Iteration %d: %d buses, %2.fs\n", i, length(routes),
                           sumIndividualTravelTimes(data, schoolID, routes))
    end
    # gc()
    GC.gc()
    return routes
 end

a = greedyCombinedIterated(data, schoolID, 1700.0, 2700.0,10, 500, 100.0; verbose = false, loglevel = 1,seconds = 10, threads = 8, allowableGap = 0.5 )
"""
    Add feasible route to a feasible set
"""
function addRoute!(routes::FeasibleRouteSet, newRoute::FeasibleRoute)
    if ! (newRoute.stopIds in routes.set)
        push!(routes.set, newRoute.stopIds)
        push!(routes.list, newRoute)
        newRouteId = length(routes.list)
        for stopId in newRoute.stopIds
            push!(routes.atStop[stopId], newRouteId)
        end
    end
end
# r = all_routes[2]
#
# newRoute = a_feasible_route
#
# routes = a_feasible_set
# routes.atStop[20:60]
#
# routeList = generateRoutes(data, schoolID, N, maxRouteTimeLower, maxRouteTimeUpper)
# routes_example = FeasibleRouteSet(data, schoolID)
# routes.atStop[20:60]



"""
    Add list of feasible routes to a feasible set
"""
function addRoute!(routes::FeasibleRouteSet, newRoutes::Vector{FeasibleRoute})
    for r in newRoutes
        addRoute!(routes, r)
    end
end


"""
    Generate N sets of Routes using the greedy heuristic.
"""
#
# maxRouteTimeLower = 2400.0
# maxRouteTimeUpper = 3600.0
# N = 3
function generateRoutes(data::SchoolBusData,
                        schoolID::Int,
                        N::Int,
                        maxRouteTimeLower::Float64=Inf,
                        maxRouteTimeUpper::Float64=Inf)
    routes = FeasibleRoute[]
    for i in 1:N
        if maxRouteTimeLower < Inf
            maxTime = (maxRouteTimeUpper - maxRouteTimeLower) * rand() + maxRouteTimeLower
        else
            maxTime = Inf
        end
        singleRoutes = greedy(data, schoolID, maxTime)
        append!(routes, [FeasibleRoute(data, schoolID,route) for route in singleRoutes])
    end
    return routes
end

# feasible_set = generateRoutes(data, 1, 3, maxRouteTimeLower, maxRouteTimeUpper )


"""
    Solves the routing problem given a set of routes. (MIP)
"""

routeList = generateRoutes(data, schoolID, 3, 1700.0, 2700.0)
routes = FeasibleRouteSet(data, schoolID)
addRoute!(routes, routeList)
function bestRoutes(data::SchoolBusData, schoolID::Int, routes::FeasibleRouteSet,
                    λ::Float64;
                    # , env::Gurobi.Env;
                     args...
                    )
    # model = Model(solver=GurobiSolver(env, Threads=getthreads(); args...))

    # model =  Model(with_optimizer(Gurobi.Optimizer; args...))
    model =  Model(with_optimizer(Cbc.Optimizer; args...))
    # The binaries, whether we choose the routes
    @variable(model, r[k in 1:length(routes.list)], Bin)
    # Minimize the number of buses first, then the travel time
    @objective(model, Min, sum(r[k] * (λ + routes.list[k].cost) for k in 1:length(routes.list)))
    # At least one route per stop
    @constraint(model, stopServed[i in 1:length(data.stops[schoolID])],
        sum(r[k] for k in routes.atStop[i]) >= 1)
    # solve(model)
     JuMP.optimize!(model)
    # selectedRoutes = [k for k in 1:length(routes.list) if getvalue(r[k]) >= 0.5]
    selectedRoutes = [k for k in 1:length(routes.list) if JuMP.value(r[k]) == 1]
    return selectedRoutes
end


selectedRoutes = bestRoutes(data, 1, routes, 100.0; seconds = 1)
"""
    Given a covering list of FeasibleRoutes, create correct route object
"""
function buildSolution(data::SchoolBusData, schoolID::Int,
                       routeSet::FeasibleRouteSet, selectedRoutes::Vector{Int})
    # retrieve the route stopid and cost for the covering model route
    routes = copy(routeSet.list[selectedRoutes])
    routesAtStop = Vector{Int}[Int[] for s in data.stops[schoolID]]
    for (routeId,r) in enumerate(routes)
        for stopId in r.stopIds
            push!(routesAtStop[stopId], routeId)
        end
    end
    # which trip covered the stop[1:END]
    for (stopId,intersectingRoutes) in enumerate(routesAtStop)
        if length(intersectingRoutes) > 1
            newRoutes = splitRoutes(data, schoolID, routes[intersectingRoutes], stopId)
            for (i, routeId) in enumerate(intersectingRoutes)
                routes[routeId] = newRoutes[i]
            end
        end
    end
    # remove the feasieble route with 0 stop because the previous deleting stop
    return [Route(i, fr.stopIds) for (i, fr) in enumerate(routes) if length(fr.stopIds) > 0]
end
#
# routeSet = routes
# routeSet.list[21]
# routes[21]
"""
    When several routes intersect in stopId, choose the best one to serve the stopId
    and remove the others
"""
function splitRoutes(data::SchoolBusData, schoolID::Int,
                     routes::Vector{FeasibleRoute}, stopId::Int)
    costs = collect(deletionCost(data, schoolID, route, stopId) for route in routes)
    # selectedRoute = indmax(costs)
    selectedRoute = argmax(costs) # max(-20, -1), lower means by reomoveing the stop, it saves more time -20
    for (routeId, route) in enumerate(routes)
        if selectedRoute != routeId
            newRouteIds = collect(s for s in route.stopIds if s != stopId)
            routes[routeId] = FeasibleRoute(newRouteIds, route.cost + costs[routeId])
        end
    end
    return routes
end
#
# stopId = 1
# intersectingRoutes = routesAtStop[stopId]
# routes =  routeSet[intersectingRoutes]
# # route = routes[intersectingRoutes][1]
#
# newRoutes =  splitRoutes(data, 1, routes, stopId)

"""
    Cost of removing a stop from a route (usually negative)
    The cost is just the difference between the old travel time and the newStop
"""

function deletionCost(data::SchoolBusData, schoolID::Int,
                      route::FeasibleRoute, stopId::Int)
    # stop = findfirst(route.stopIds, stopId)
    stop = findfirst(isequal(stopId), route.stopIds)
    length(route.stopIds) <= 1 && return (-Inf)
    stops = data.stops[schoolID]
    school = data.schools[schoolID]
    # remove travel time of that particular stop
    cost = 0.
    for id = stop:(length(route.stopIds)-1)
        # cost -= traveltime(data, stops[route.stopIds[id]], stops[route.stopIds[id+1]])
        # cost -= stopTime(data, stops[route.stopIds[id+1]])
        global cost -= traveltime(data, stops[route.stopIds[id]], stops[route.stopIds[id+1]])
        global cost -= stopTime(data, stops[route.stopIds[id+1]])
    end
    cost -= traveltime(data, stops[route.stopIds[end]], school)
    # remove travel time effect on other stops
    if stop == length(route.stopIds)
        global cost -= (traveltime(data, stops[route.stopIds[stop-1]], stops[stopId]) +
                 traveltime(data, stops[stopId], school) +
                 stopTime(data, stops[stopId]) -
                 traveltime(data, stops[route.stopIds[stop-1]], school)) *
                (length(route.stopIds) - 1)
    elseif stop > 1
        global cost -= (traveltime(data, stops[route.stopIds[stop-1]], stops[stopId]) +
                 traveltime(data, stops[stopId], stops[route.stopIds[stop+1]]) +
                 stopTime(data, stops[stopId]) -
                 traveltime(data, stops[route.stopIds[stop-1]], stops[route.stopIds[stop+1]]))*
                (stop - 1)
    end
    return cost
end

# deletionCost(data, 1, routes[1], 1 )


"""
    Contains the parameters that were used to compute a particular scenario
"""
struct ScenarioParameters
    "Maximum time of routes - lower end of interval"
    maxRouteTimeLower::Float64
    "Maximum time of routes - upper end of interval"
    maxRouteTimeUpper::Float64
    "Number of greedy solutions optimized over"
    nGreedy::Int
    "Tradeoff parameter in optimization"
    λ::Float64
    "Number of iterations for greedy combined iterated"
    nIterations::Int
end
"""
    Constructor for ScenarioParameters, with keyword parameters for readability
"""
function ScenarioParameters(;maxRouteTimeLower=Inf,
                            maxRouteTimeUpper=Inf,
                            nGreedy::Int=10,
                            λ=5e3,
                            nIterations::Int=10)
    return ScenarioParameters(maxRouteTimeLower, maxRouteTimeUpper, nGreedy, λ, nIterations)
end

λs = [1e2, 5e2, 1e3, 2e3, 5e3];

maxtime = 2700.0
nGreedy  = 10
params = [ScenarioParameters(maxRouteTimeLower=maxtime-2000,
                                 maxRouteTimeUpper=maxtime+2000,
                                 nGreedy=nGreedy, λ=λ,
                                 nIterations=10) for λ in λs]
scenariolist = computescenarios(data, params;
                            ## Gurobi arguments
                                 # OutputFlag=0,
                                 # MIPGap=ifelse(maxtime < 4000, 1e-4, 0.05),
                                 # TimeLimit=ifelse(maxtime < 4000, 90, 30)
                            ## Cbc arguments
                                 logLevel=1,
                                 allowableGap=ifelse(maxtime < 4000, 1e-4, 0.05),
                                 seconds=ifelse(maxtime < 4000, 90, 30),
                                 threads = 8
                                 );
"""
    Compute scenario
"""
function getscenario(data, scenarioinfo; args...)
    school, scenarioid, params = scenarioinfo
    routes = greedyCombinedIterated(data, school.id,
                                    params.maxRouteTimeLower, params.maxRouteTimeUpper,
                                    params.nGreedy, params.nIterations, params.λ;
                                    verbose=false, args...)
    return Scenario(school.id, scenarioid, collect(eachindex(routes))), routes
end

"""
    Compute multiple scenarios
    tocompute:
30-element Array{Tuple{School,Int64,ScenarioParameters},1}:
 (School 5 - 200005, 4, ScenarioParameters(700.0, 4700.0, 10, 2000.0, 10))
 (School 4 - 200004, 3, ScenarioParameters(700.0, 4700.0, 10, 1000.0, 10))
 (School 4 - 200004, 2, ScenarioParameters(700.0, 4700.0, 10, 500.0, 10))
 (School 5 - 200005, 2, ScenarioParameters(700.0, 4700.0, 10, 500.0, 10))
"""
function computescenarios(data, params; args...)
    tocompute = shuffle!(vec([(school,paramid,param) for school in data.schools, (paramid,param) in enumerate(params)]))
    results = Tuple{Scenario,Vector{Route}}[]
    @showprogress for scenarioinfo in tocompute
        push!(results, getscenario(data, scenarioinfo; args...))
    end
    return results
end
function computescenariosparallel(data, params; args...)
    tocompute = shuffle!(vec([(school,paramid,param) for school in data.schools, (paramid,param) in enumerate(params)]))
    results = Tuple{Scenario,Vector{Route}}[]
    results = pmap(x->getscenario(data,x;args...), tocompute)
    return results
end

"""
    Put the scenarios together
"""
function loadroutingscenarios!(data, scenariolist)
    scenarios = [Scenario[] for school in data.schools]
    routes = [Route[] for school in data.schools] # each school has a set of routes
    ids = vec([(scenario.school, scenario.id) for (scenario, routelist) in scenariolist])

    for k in sortperm(ids)
        (scenario, routelist) = scenariolist[k]
        for (i,route) in enumerate(routelist)
            idx = findfirst(x -> x.stops == route.stops, routes[scenario.school])
            # if idx == 0 # need to add new route and update its id
            if idx == nothing # need to add new route and update its id
                scenario.routeIDs[i] = length(routes[scenario.school]) + 1
                push!(routes[scenario.school],
                      Route(length(routes[scenario.school]) + 1, route.stops))
            else # need to update reference in the scenario
                scenario.routeIDs[i] = idx
            end
        end
        push!(scenarios[scenario.school], scenario)
    end
    data.scenarios = scenarios
    data.routes = routes
    data.withRoutingScenarios = true
    return data
end

k = 11
i = 1
route = Route(1, [67])
data = loadroutingscenarios!(data, scenariolist);

data.scenarios[1]
data.routes[1][9]
