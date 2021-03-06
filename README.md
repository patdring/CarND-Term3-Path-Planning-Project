# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

### Goals
In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.

[//]: # (Image References)
[image1]: ./images/no_car_in_ahead_2.png
[image2]: ./images/no_car_in_ahead.png
[image3]: ./images/overtaking_red_car.png
[image4]: ./images/overtook_red_car.png
[image5]: ./images/path_planner.png
[image6]: ./images/target_dist.png

#### Requirements summary for valid trajectories
* The car is able to drive at least 4.32 miles without incident

![alt text][image2]

* The car drives according to the speed limit.

![alt text][image1]

* Max Acceleration and Jerk are not Exceeded
* Car does not have collisions
* The car stays in its lane, except for the time between changing lanes

![alt text][image3]

* The car is able to change lanes

![alt text][image4]

### Model Documentation

#### Overview
Path Planning includes the six components Behaviour, Prediction, Localization, Trajectory, Sensor Fusion and Motion Control. The components marked green in the following diagram are to be implemented in this project.

![alt text][image5]

#### PREDICTION
This component tries to predict where other vehicles might be in the near future.

In addition to predicting at which position S a vehicle will be (a possible trajectory), a possible change of lane is also predicted. For this purpose a Gaussian Naive Bayes classifier is used. The corresponding exercise from the "Prediction" module was used as code base and training data. So a hybrid approach is used that combines a model-based and a data driven approach.

Example data set of a vehicle which is in 30m range, behind the ego car and in lane 2 so right behind us and it seems to be that it's keeping lane.

```
id    5
x     1082.08
y     1170.37
vx    19.5193
vy    7.4421
s     301.402
s_rel -25.1249
d     10.149
s.    -0.0990217
d.    0.0502717
cl    2
pl    keep
---------------------------------------
No. of other vehicles in range of +/- 30m: 1
No. of (OVs) lane changes detected:  0
is_vehicle_left  0
is_vehicle_right 1
is_vehicle_front 0
``` 

```
- `s_rel`: Relative distance to ego car. Negative means it's behind our ego car
- `s.`: s dot, change of S to time t. Neccessary as param for Gaussian Naive Bayes classifier
- `d.`: d dot, change of d to time t. Neccessary as param for Gaussian Naive Bayes classifier 
- `cl`: Current lane 
- `pl`: Predicted lane in future for this vehicle returnd by Gaussian Naive Bayes classifier
```

Finally a vector containing all vehicles will be reduced to the vehicles within a range of +/- 30m. Flags are used to record where these vehicles are or will be located relative to your own vehicle.

#### BEHAVIOUR
This component determines where the ego vehicle must be in the future. Specifically, the maneuvers "Keep Lane", "Change Left" or "Change Right". In general, of course, more like stopping at traffic lights or intersections, roundabouts, emergency stops and more.

Cost functions are used for complex planning. Cost functions are used to select the best trajectory but for the mentioned simple manoeuvres this is omitted and the already mentioned flags are used.

The following code snippet now implements the behavior of the ego vehicle.

```
if (is_vehicle_front) { 
  // case: no vehicle left and it's possible to change
  if (!is_vehicle_left && lane > 0) {          
    lane--; 
  // case: no vehicle right and it's possible to change
  } else if (!is_vehicle_right && lane != 2) {
    lane++; 
  // case: it's not possible to change lanes so reduce speed
  } else {
    speed_diff -= 0.224;
  }
} else {
  // case: no vehicle in front so increase speed to the allowed limit
  if (ref_vel < 49.5) {
    speed_diff += 0.224;
  }
}
```

#### TRAJECTORY
The creation of a trajectory is based on the idea presented in the course video "Project Q&A".

1. Save size of points of a previous path
2. Creating reference values of x,y and yaw
3. If there are not enough points of a previous path, we use the data of the ego car to find a previous point and store these two in a vector
4. If points of a previous path are existing then the last two points are simply added to the vectors and the last point is set as reference
5. Now 3 points are defined which are in the future. These points are determined with the waypoints of the map
6. Add remaining previous points to the vectors
7. A spline is used for the trajectory generation. The library for this is easy to use and there are no dependencies
    - https://en.wikipedia.org/wiki/Spline_(mathematics)
    - http://kluge.in-chemnitz.de/opensource/spline/
7. Find the all spline points next 30m so that spacing the way that ego car can travel at desired speed
   ```
   double target_x = 30.0;
   double target_y = s(target_x);
   double target_dist = sqrt(target_x*target_x + target_y*target_y);
   ```
   ![alt text][image6]
8. Calculate the spline points (50-previous_points) from start to horizon y points
9. Transfer of the final control values to the simulator

### Simulator.
You can download the Term3 Simulator which contains the Path Planning Project from the [releases tab (https://github.com/udacity/self-driving-car-sim/releases/tag/T3_v1.2).

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

## Tips

A really helpful resource for doing this project and creating smooth trajectories was using http://kluge.in-chemnitz.de/opensource/spline/, the spline function is in a single hearder file is really easy to use.

---

## Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!


## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./

## How to write a README
A well written README file can enhance your project and portfolio.  Develop your abilities to create professional README files by completing [this free course](https://www.udacity.com/course/writing-readmes--ud777).

