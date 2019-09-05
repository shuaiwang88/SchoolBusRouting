# Data Requirement

## Default parameters

| Name                  | Default Value |
|-----------------------|:-------------:|
| bus_capacity          |       66      |
| max_time_on_bus       |      3600     |
| constant_stop_time    |       30      |
| stop_time_per_student |       5       |
| velocity              |      20       |
| metric                |   MANHATTAN   |


## Point (lat/long)

| Name | Type  |
|------|-------|
| x    | Float |
| y    | Float |


## School 

| Name       | Type  | Comment                                                                  |
|------------|-------|--------------------------------------------------------------------------|
| School Id  | Int   | index of the school                                                      |
| OriginalId | Int   | Unique school identifier                                                 |
| Postion    | Point | position(x,y)                                                            |
| dwelltime  | Float | dwell time/drop-off time                                                 |
| starttime  | Float | School start time(in sec from 12am): deadline for bus to finish drop-off |
| start_tw   | Float | time window start                                                        |
| end_tw     | Float | time window end                                                          |
