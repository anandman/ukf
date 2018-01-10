# Unscented Kalman Filter

This code implements an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy LIDAR and RADAR measurements. 

---

## Dependencies
* [Udacity Self-Driving Car Simulator](https://github.com/udacity/self-driving-car-sim/releases)
* [uWebSocketIO](https://github.com/uWebSockets/uWebSockets)
    * Linux: run the [`install-ubuntu.sh`](install-ubuntu.sh) script in the repository
    * Mac: run the [`install-mac.sh`](install-mac.sh) script in the repository
    * Windows: use either Docker, VMware, or even [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
* cmake >= 3.5
    * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
    * Linux: make is installed by default on most Linux distros
    * Mac: install [Xcode command line tools](https://developer.apple.com/xcode/features/) to get make (type `xcode-select --install` in Terminal)
    * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
    * Linux: gcc / g++ is installed by default on most Linux distros
    * Mac: same deal as make - install [Xcode command line tools](https://developer.apple.com/xcode/features/)
    * Windows: recommend using [MinGW](http://www.mingw.org/)

## Build and Run Instructions

1. Clone this repo
2. Install uWebSocketIO as indicated [above](#dependencies)
3. Make a build directory: `mkdir build && cd build`
4. Compile: `cmake .. && make`
5. Run it: `./UnscentedKF`
6. Start Udacity simulator application
    1. Select resolution, graphics quality, etc., and press "Play!" 
    2. Select "Project 1/2: EKF and UKF"
    3. Pick Dataset and press "Start"

## uWebSocketIO Protocol Usage

If you wish to modify the communication between the simulator and filter, here is the main protocol that main.cpp uses for uWebSocketIO.

INPUT: values provided by the simulator to the C++ program

["sensor_measurement"] => the measurment that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

