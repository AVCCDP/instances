# AVCCDP
This repository includes the instances and source code used in our paper "Meals on Autonomous Wheels: Delivery with Vehicle-Customer Coordination".

With the rapid advancement of autonomous driving technology, autonomous meal delivery has emerged as a promising solution for efficient urban logistics. This paper studies the Autonomous Vehicle-Customer Coordination Delivery Problem (AVCCDP). Specifically, we introduce a vehicle-customer coordination mechanism in which autonomous vehicles serve as mobile lockers (secure, temporary pickup points). This mechanism enables flexible parking and guided customer movement to improve delivery efficiency while maintaining user convenience. Given the time-sensitive nature of meal delivery and the need to encourage customer movement, it is essential to jointly optimize delivery routes and customer movements. We formalize this challenge as a nondeterministic polynomial-time (NP)-hard problem. Furthermore, the AVCCDP provides a generalized framework that is applicable to both the order batching policy (OBP) and the order postponement policy (OPP).

## Instances
The [repository](./repository) folder contains all data instances used in this research, which are divided into two main categories: offline instances and online instances.

### (i) Offline Instances
- All offline networks are modified based on the Solomon benchmark for VRPTW. The letters C, R, and RC refer to clustered distribution networks, random distribution networks, and clustered-random mixed distribution networks, respectively. The number of customers and the case number are also included in the filename. For example, "C1_N20" indicates that this instance is the first clustered distribution case with 20 customer points. The instance file format is described as follows:
  - StringID: Node identifier (D represents the depot; C represents the customer).
  - x, y: Node coordinates.
  - profit: The maximum profit obtainable at this node.
  - demand: Node demand.
  - delayRatio: Linear decay rate of profit with arrival time delay.
  - Tmax: The latest time limit for the vehicle to return to the depot.

### (ii) Online Instances
- To simulate real-world delivery scenarios, we generated 5 sets of online instances based on the order flow of recommended merchants on the Uber Eats homepage during peak hours. The order arrivals follow a Poisson process. Each set of instances covers one peak hour and is divided into 12 scheduling periods (5 minutes per period). To comprehensively examine the effectiveness of the coordination mechanism, we designed the following comparison dimensions under a rolling time horizon framework:
  - Scheduling Strategy
    - OBP (Order Batching Policy): All orders in the current period must be dispatched, even if they may exceed the time limit.
    - OPP (Order Postponement Policy): This policy allows orders that are expected to exceed the time limit to be postponed to subsequent periods.
  - Load Scenario
    - Low Load Rate (LLR): At least 2 vehicles are assigned per period, with a maximum of 5 orders per vehicle.
    - High Load Rate (HLR): Only 1 vehicle is assigned per period, with a maximum of 10 orders per vehicle.
    - In both scenarios, additional vehicles are dispatched only when the number of orders exceeds the capacity of a single vehicle.

- The filename of each instance indicates the case number and scheduling period. For example, "Case1_T1" represents the orders that need to be assigned in the first scheduling period of the first case. The following file suffixes indicate the corresponding matrix data, which was obtained via the Google Distance Matrix API:
  - _drivingDis: Driving distance matrix between nodes.
  - _drivingTime: Driving time matrix between nodes.
  - _walkingDis: Walking distance matrix between nodes.
  - _walkingTime: Walking time matrix between nodes.

- To distinguish order information under different scheduling strategies, the filenames may include the following suffixes:
  - _origin: This suffix indicates the original order data under the OBP strategy during this period.
  - _move: This suffix indicates the order data under the OPP strategy with the coordination mechanism enabled during this period.
 
- Regarding the file format, online instances follow the same structure as offline instances, with one addition: the "StringID" field includes a new node identifier "S" to represent public parking spots.

## Algorithm
We developed a branch-and-price framework equipped with a customized re-labeling algorithm. The corresponding code scripts are provided in the [src](./repository/src) folder.
