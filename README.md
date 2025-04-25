# Demand-responsive Scheduler

This project provides functionalities to create and solve simple real-world demand-responsive problem instances.
Specifically, it defines instances of a dynamic pickup and delivery problem with time windows.

The implementation of this project is based on the work of 
[Mark E.T. Horn (2002)](https://www.sciencedirect.com/science/article/pii/S0968090X01000031) and thoroughly described in
our own article [Martí et al. (2024)](https://doi.org/10.2298/CSIS230115074M)

### Demand-responsive transportation.
Fully-dynamic demand-responsive transportation services solve transportation demands by dynamically creating routes for
their fleet vehicles. Transportation demand is represented by customer requests. Each request defines a trip of a number
of passengers between an origin and a destination stop and that shall be completed within a certain time window. 
An example of a request could be "Pickup 2 customers at stop A after 8:30 and dropoff at stop B no later than 9:30".

A problem instance is solved by **inserting customer trips** (represented in the requests) **into 
vehicle itineraries** (represented by fleet vehicles). Each vehicle will follow its own itinerary:
the list of stops it will visit during its operation, ordered in time. Inserting a trip into an itinerary implies 
planning a **visit to the trip's origin stop** and a posterior **visit to the trip's destination stop**. Both visits 
must **preserve the temporal restrictions** defined by the time window of the customer's request. In addition, 
**the capacity** of the vehicle serving the request **must not be exceeded** at any point of its operations.  



## Scheduling algorithm
A fleet of vehicles provides displacement services with a variable capacity. Each vehicle will follow its own itinerary: 
the list of stops it will visit during its operation, ordered in time. We assume that users of the system issue travel 
requests through an application. A travel request indicates the location and time window in a simple manner, 
such as ``Pickup at stop A after 8:30, and dropoff at B by 9:00".
Our system is managed by a centralized scheduler which allocates each travel request to a vehicle's itinerary.
The scheduler allocates the requests to itineraries such that the system-wide objective function is optimized. 
## Input data
The project requires several input files to create and solve demand-responsive transportation problem instances. 
Input data is composed by the **Stops file**, the **Routes file** and a specific **Problem Configuration file**.

> Note: The Stops file must be **manually created by the user**. Then, the project provides functionalities to generate the 
> Routes file according to a given Stops file with the `routeCalculator` module. Having the Stops and the Routes file, 
> different problem instances can be generated through the `vehicleGenerator` and `customerGenerator` modules.

#### Stops file
- Stops file (.json). Real-world locations that represent the stops of the transportation service. 
The format of the Stops file follows that of [GeoJSON](https://geojson.org/) files. Each stop is defined by:
  - `id`
  - `coordinates`
  - (optional) `properties`
- Example of a Stops file:
```
{
    "type": "FeatureCollection",
    "features": [
        {
            "id": "0",
            "type": "Feature",
            "properties": {},
            "geometry": {
                "type": "Point",
                "coordinates": [
                    -0.3005497742,
                    39.1969752628
                ]
            }
        },
        ...
    ]
}
```
#### Routes file
- Routes file (.json). Dictionary indexed by pairs of stop coordinates that defines the route between each pair 
of different stops. The format of the Routes file follows that of the replied given by the routing server, which
is implemented with [Open Source Routing Machine (OSRM)](https://project-osrm.org/). Each route is defined by:
  - `path` (list of coordinates of the driving route that connects the stops)
  - `distance` (of the route that connects the stops, in meters)
  - `duration` (estimation of time taken by a transport to traverse the route, in seconds)
- Example of a Routes file:
```
{
    "(-0.3005497742, 39.1969752628):(-0.305559128, 39.2001829203)": {
        "path": [ [39.196931, -0.300586], ..., [-0.305559128, 39.2001829203] ], 
        "distance": 588.6, 
        "duration": 52.7
    },
    ...
}
```
#### Problem Configuration file
- Configuration file (.json). Dictionary describing the vehicle fleet and customer demand
of an experiment. All agents must be localised in one of the stops of the Stop file. Each configuration file defines
a problem instance. The format of the configuration files is derived from the 
[SimFleet](https://github.com/javipalanca/simfleet) configuration files
  - `"transports"`: List of dictionaries, one per fleet vehicle, defining:
    - `name`,
    - `position` (origin stop), `destination` (final stop), (coordinates)
    - `capacity`,
    - `speed`,
    - `start_time` (beginning of shift at position), `end_time` (end of shift at destination)

  - `"customers"`: List of dictionaries, one per customer request, defining
    - `name`,
    - `position` (origin stop), `destination` (destination stop), (coordinates)
    - `npass` (number of passengers travelling together)
    - `issue_time` (time at which the customer issues their request)
    - `origin_time_ini`, `origin_time_end` (temporal window for pickup at origin stop)
    - `destination_time_ini`, `destination_time_end` (temporal window for arrival at the destination)
- Example of a Configuration file:
```
{
    "transports": [
        {
            "name": "veh_01"
            "position": [
                39.4464599991,
                -0.326572
            ],
            "destination": [
                39.4464599991,
                -0.326572
            ],
            "capacity": 8,
            "speed": 70,
            "start_time": 0,
            "end_time": 960,            
        },
        ...
    ],
    "customers": [
        {
            "name": "request_001"
            "position": [
                39.3624981331,
                -0.4577106041
            ],
            "destination": [
                39.349443,
                -0.323432
            ],
            "npass": 1,
            "issue_time": 61,
            "origin_time_ini": 66,
            "origin_time_end": 86,
            "destination_time_ini": 88.13,
            "destination_time_end": 125.27,           
        },
        ...
    ]
}    
```


## Problem instance generation
Problem instances are generated departing from a 

### Route calculation
`routeCalculator.py`

### Vehicle generator
`vehicleGenerator.py`

### Customer generator
`customerGenerator.py`

### Run experiments

## Solving a problem instance