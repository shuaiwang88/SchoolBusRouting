###################################################
## load.jl
##     Loads benchmarks from Park, Tae and Kim (2011)
## Authors: Arthur Delarue, Sébastien Martin, 2018
###################################################
# include("problem.jl")
using Random
using CSV
using DataFrames
"""
    Convert weird coordinates from synthetic benchmarks to lat and lon
        Assumes a random center
"""
function parseCoordinates(x::Real, y::Real; center="Default")
    # for these benchmarks we'll keep the weird units
    return Point(x, y)
end

parseCoordinates(111, 222.33)

"""
    Helper function to parse time from synthetic text files
"""
function parseTime(time)
    s = string(time)
    minutes = parse(Int8,s[end-1:end]) * 60
    hours = parse(Int8,s[1:end-2]) * 3600
    return hours+minutes
end

parseTime("0630")
"""
    Load a tab-separated file containing the school data
"""
function loadSchoolsReduced(schoolsFileName::AbstractString,
                            maxEffect::Real=Inf,
                            randomStart::Bool=false, seed::Int=-1,
                            spreadStart::Bool=false)
    if seed < 0
        # srand()
        Random.seed!()
    else
        #srand(seed)
        Random.seed!(seed)
    end
    schoolData  = CSV.read(schoolsFileName)
    #= schoolData  = CSV.read(schoolsFileName, delim="\t", DataFrame) =#

    if spreadStart
        arrivalTimes = spreadBellTimes(schoolData, maxEffect)
    end
    schools = School[]
    for i in 1:nrow(schoolData)
        id = length(schools) + 1
        #= originalId = get(schoolData[i, :ID]) =#
        originalId = (schoolData[i, :ID])
        #= position = parseCoordinates(get(schoolData[i, :X]), get(schoolData[i, :Y])) =#
        position = parseCoordinates((schoolData[i, :X]), (schoolData[i, :Y]))
        dwelltime = 154.4
        #= intervalstart = parsetime(get(schooldata[i, :amearly])) =#
        #= intervalend = parsetime(get(schooldata[i, :amlate])) =#
        intervalStart = parseTime((schoolData[i, :AMEARLY]))
        intervalEnd = parseTime((schoolData[i, :AMLATE]))
        if randomStart # randomly select start time in allowed window
            starttime = intervalStart + (intervalEnd-intervalStart)*rand()
        elseif spreadStart
            starttime = arrivalTimes[i]
        else
            starttime = intervalStart
        end
        push!(schools, School(id, originalId, position, dwelltime, starttime,
                              intervalStart, intervalEnd))
    end
    return schools
end
#
# test_school = loadSchoolsReduced("../data/input/CSCB01/Schools.txt")
# test_school[1].position
#
function spreadBellTimes(schoolData::DataFrame, maxEffect::Real)
    intervalStart = [parseTime(get(schoolData[i,:AMEARLY])) for i=1:nrow(schoolData)]
    intervalEnd = [parseTime(get(schoolData[i,:AMLATE])) for i=1:nrow(schoolData)]
    intervalStart = 7.5 * 3600 + 2 * (intervalStart - minimum(intervalStart))/
    							 (maximum(intervalStart) - minimum(intervalStart))
    return intervalStart
    # model = Model(solver=GurobiSolver(OutputFlag=0))
    # @variable(model, intervalStart[i] <= belltime[i=1:nrow(schoolData)] <= intervalEnd[i])
    # # distances
    # @variable(model, 0 <= d[i=1:nrow(schoolData), j=1:nrow(schoolData)] <=
    #                  maximum(intervalEnd)-minimum(intervalStart))
    # @constraint(model, [i=1:nrow(schoolData), j=1:nrow(schoolData);
    #                     0 < intervalStart[j] - intervalEnd[i] < maxEffect],
    #             d[i,j] <= belltime[j] - belltime[i])
    # @objective(model, Max, sum(d))
    # solve(model)
    # return getvalue(belltime)
end

"""
    Create a yard in the center of the district, with 200 full buses
"""
function syntheticYards()
	return [Yard(1, parseCoordinates(105_600., 105_600.))]
end

"""
    Load bus stops with students
"""
function loadPreComputedStops(stopsFileName::AbstractString, schools::Vector{School})
    #= stopData = CSV.read(stopsFileName, delim="\t", DataFrame) =#
    stopData = CSV.read(stopsFileName)
    schoolIdMap = Dict(school.originalId => school.id for school in schools)
    stops = [Stop[] for school in schools]
    for i = 1:nrow(stopData)
        originalId = (stopData[i, :ID])
        schoolId = schoolIdMap[(stopData[i,:EP_ID])]
        position = parseCoordinates((stopData[i,:X_COORD]), (stopData[i,:Y_COORD]))
        nStudents = (stopData[i,:STUDENT_COUNT])
        push!(stops[schoolId],
              Stop(length(stops[schoolId])+1, i, originalId, schoolId, position, nStudents))
    end
    return stops
end
# 
# bus_stops = loadPreComputedStops("../data/input/CSCB01/Stops.txt", test_school)
# #= stopData = CSV.read("../data/input/CSCB01/Stops.txt") =#
#
# test_school
#

"""
    Load synthetic benchmark dataset
"""
function loadSyntheticBenchmark(schoolsFile::AbstractString, stopsFile::AbstractString;
                                randomStart::Bool=false, seed::Int=-1,
                                spreadStart::Bool=false, maxEffect::Real=Inf)
	data = SchoolBusData()
    data.schools = loadSchoolsReduced(schoolsFile, maxEffect, randomStart, seed, spreadStart)
    data.yards = syntheticYards()
    data.stops = loadPreComputedStops(stopsFile, data.schools)
    # update parameters
    data.params.bus_capacity                         = 66
    data.params.max_time_on_bus                      = 2700.
    data.params.constant_stop_time                   = 19.
    data.params.stop_time_per_student                = 2.6
    data.params.velocity                             = 29.3333333#20 * 0.44704
    data.params.metric                               = MANHATTAN
    data.withBaseData = true
    return data
end
