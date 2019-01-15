#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "classifier.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

vector<vector<double> > Load_State(string file_name)
{
    ifstream in_state_(file_name.c_str(), ifstream::in);
    vector< vector<double >> state_out;
    string line;
    
    
    while (getline(in_state_, line)) 
    {
        istringstream iss(line);
    	vector<double> x_coord;
    	
    	string token;
    	while( getline(iss,token,','))
    	{
    	    x_coord.push_back(stod(token));
    	}
    	state_out.push_back(x_coord);
    }
    return state_out;
}
vector<string> Load_Label(string file_name)
{
    ifstream in_label_(file_name.c_str(), ifstream::in);
    vector< string > label_out;
    string line;
    while (getline(in_label_, line)) 
    {
    	istringstream iss(line);
    	string label;
	    iss >> label;
    
	    label_out.push_back(label);
    }
    return label_out;
    
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// start in lane 1
int lane = 1;
// reference velocity to target [mph]
double ref_vel = 0.0; 
int total_detected_lane_changes = 0;

struct vehicle {
  int id;
  double x;
  double y; 
  double vx; 
  double vy;
  double speed;
  double s;
  double s_rel;
  double pred_s;
  double d;
  double s_dot;
  double d_dot;
  string pred_lane;
  int curr_lane;
};
  
vector<vehicle> prev_other_vehicles;

GNB gnb = GNB();

std::chrono::steady_clock::time_point newT = std::chrono::steady_clock::now();
std::chrono::steady_clock::time_point oldT = newT;

/*
typedef enum
{
  KEEP_LANE = 0,
  CHANGE_LEFT,
  CHANGE_RIGHT,
  PREPARE_CHANGE_LEFT,
  PREPARE_CHANGE_RIGHT  
} path_planner_state_t;

path_planner_state_t pp_state = KEEP_LANE;
*/

int main() {
  uWS::Hub h;
 
  vector< vector<double> > X_train = Load_State("../data/train_states.txt");   
  vector< string > Y_train  = Load_Label("../data/train_labels.txt");

  cout << "Start Naive Bayes Classifier training ..." << endl;
  cout << "X_train number of elements " << X_train.size() << endl;
  cout << "X_train element size " << X_train[0].size() << endl;
  cout << "Y_train number of elements " << Y_train.size() << endl;

  gnb.train(X_train, Y_train);
  cout << "finished!" << endl;
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  
  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      newT = std::chrono::steady_clock::now();
      auto time_diff = std::chrono::duration_cast<std::chrono::milliseconds>(newT - oldT).count();
      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
                 
          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];
          newT = std::chrono::steady_clock::now();
          
          int prev_size = previous_path_x.size();
          if(prev_size > 0) {
            car_s = end_path_s;
          }
          
          json msgJson;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          vector<vehicle> other_vehicles;
                  
          for (int i=0; i < sensor_fusion.size(); i++) {
            vehicle v;
            
            v.id = sensor_fusion[i][0];
            v.x = sensor_fusion[i][1];
            v.y = sensor_fusion[i][2];
            v.vx = sensor_fusion[i][3]; 
            v.vy = sensor_fusion[i][4];
            v.speed = sqrt(v.vx*v.vx+v.vy*v.vy);
            v.s = sensor_fusion[i][5];
            v.s_rel = v.s - car_s;
            v.d = sensor_fusion[i][6];
            v.s_dot = 0;
            v.d_dot = 0;
            v.pred_lane = "";
            v.pred_s = v.s + ((double)prev_size*0.02*v.speed);
                                  
            if (v.d > 0 && v.d < 4) {
              v.curr_lane = 0;
            } else if (v.d > 4 && v.d < 8) {
              v.curr_lane = 1;
            } else if (v.d > 8 and v.d < 12) {
              v.curr_lane = 2;
            }
                
            if (prev_other_vehicles.size() != 0) {
              
              v.s_dot = (v.s_rel - prev_other_vehicles[i].s_rel) / ((double)prev_size*0.02); 
              v.d_dot = (v.d - prev_other_vehicles[i].d) / ((double)prev_size*0.02);
                    
              // observation is a tuple with 4 values: s, d, s_dot and d_dot.
              if (abs(v.s_dot) < 13.0 && ( v.d_dot > -2.5 && v.d_dot < 2.5)) {
                vector<double> coords;
               
                coords.push_back(abs(v.s_rel));
                coords.push_back(v.d);
                coords.push_back(abs(v.s_dot));
                coords.push_back(v.d_dot);
                
                v.pred_lane = gnb.predict(coords);
                
                if (v.pred_lane == "left" && lane > 0) {
                   v.curr_lane--;
                }               
                if (v.pred_lane == "right" && lane != 2) {
                   v.curr_lane++;
                }
                
              }
            }           
            other_vehicles.push_back(v);            
          }
          //cout << "time duration: "<< std::chrono::duration_cast<std::chrono::milliseconds>(newT - oldT).count()/100.0 << endl;
          
          prev_other_vehicles.clear();
         
          for (int i=0; i < other_vehicles.size(); i++) {
            prev_other_vehicles.push_back(other_vehicles[i]);
          }
          
          vector<vehicle> cons_other_vehicles;
                     
          for (int i=0; i < other_vehicles.size(); i++) {       
            if (other_vehicles[i].s_rel < -30.0 || other_vehicles[i].s_rel > 30.0) {
              continue;
            }
            cons_other_vehicles.push_back(other_vehicles[i]); 
          }   
            
           
          /*
          cout << "x   " << car_x << endl;
          cout << "y   " << car_y << endl;
          cout << "s   " << car_s << endl;
          cout << "d   " << car_d << endl;
          cout << "yaw " << car_yaw << endl;
          cout << "s   " << car_speed << endl;
          
          cout << "+++++++++++++++++++++++++" << endl;
          */
          
          
          for (int i=0; i < cons_other_vehicles.size(); i++) {
            cout << "id    " << cons_other_vehicles[i].id << endl;
            cout << "x     " << cons_other_vehicles[i].x << endl;
            cout << "y     " << cons_other_vehicles[i].y << endl;
            cout << "vx    " << cons_other_vehicles[i].vx << endl;
            cout << "vy    " << cons_other_vehicles[i].vy << endl;
            cout << "s     " << cons_other_vehicles[i].s << endl;
            cout << "s_rel " << cons_other_vehicles[i].s_rel << endl;
            cout << "d     " << cons_other_vehicles[i].d << endl;
            cout << "s.    " << cons_other_vehicles[i].s_dot << endl;
            cout << "d.    " << cons_other_vehicles[i].d_dot << endl;
            cout << "cl    " << cons_other_vehicles[i].curr_lane << endl;
            cout << "pl    " << cons_other_vehicles[i].pred_lane << endl;
            cout << "---------------------------------------" << endl;
            
            if (cons_other_vehicles[i].pred_lane == "right" || cons_other_vehicles[i].pred_lane == "left") {
              total_detected_lane_changes++;
            }
          }
          
          /*
          double vehicle_min_dist = 30.0;
          double target_vel = 49.5;
          
          // path planner finite state machine                 
          for (int i=0; i < cons_other_vehicles.size(); i++) { 
            
            switch (pp_state) {
              case KEEP_LANE:                                           
                // is car in same lane and in front of our ego vehicle?
                double curr_dist = cons_other_vehicles[i].pred_s - car_s;
                if (cons_other_vehicles[i].curr_lane == lane && cons_other_vehicles[i].pred_s > car_s && curr_dist < vehicle_min_dist) {               
                  // search car wich is directly in front of our ego car and takeover its speed           
                  vehicle_min_dist = curr_dist;
                  target_vel = cons_other_vehicles[i].speed;
                 
                  if ((cons_other_vehicles[i].curr_lane - lane) == -1 && lane > 0) {
                    pp_state = PREPARE_CHANGE_LEFT;              
                    break;
                  }
 
                }                                                                                  
                pp_state = KEEP_LANE;              
                break;
              case CHANGE_LEFT:
                break;
              case CHANGE_RIGHT:
                break;
              case PREPARE_CHANGE_LEFT:
                break;
              case PREPARE_CHANGE_RIGHT:
                break;
            }
          }
          
          if (ref_vel < target_vel) {
            ref_vel += 1.12;
          }
              
          if (ref_vel > target_vel) {
            ref_vel -= 1.12;
          }  
        
          cout << "No. of other vehicles in range of +/- 30m: " << cons_other_vehicles.size() << endl;
          cout << "No. of (OVs) lane changes detected:  " << total_detected_lane_changes << endl;
          */
                
          bool is_vehicle_left = 0;
          bool is_vehicle_right = 0;
          bool is_vehicle_front = 0;
                  
          for (int i=0; i < cons_other_vehicles.size(); i++) {  
            if(cons_other_vehicles[i].curr_lane == lane) {
              is_vehicle_front |= cons_other_vehicles[i].pred_s > car_s && (cons_other_vehicles[i].pred_s - car_s) < 30;           
            } else if((cons_other_vehicles[i].curr_lane - lane) == -1) {                            
              is_vehicle_left |= (car_s+30) > cons_other_vehicles[i].pred_s  && (car_s-30) < cons_other_vehicles[i].pred_s;
            } else if((cons_other_vehicles[i].curr_lane - lane) == 1) {
              is_vehicle_right |= (car_s+30) > cons_other_vehicles[i].pred_s  && (car_s-30) < cons_other_vehicles[i].pred_s;
            }        
          }

          cout << "No. of other vehicles in range of +/- 30m: " << cons_other_vehicles.size() << endl;
          cout << "No. of (OVs) lane changes detected:  " << total_detected_lane_changes << endl;
          cout << "is_vehicle_left  " << is_vehicle_left << endl;
          cout << "is_vehicle_right " << is_vehicle_right << endl;
          cout << "is_vehicle_front " << is_vehicle_front << endl;
          
          double speed_diff = 0;
         
          if (is_vehicle_front) { 
            // case no vehicle left and it's possible to change
            if (!is_vehicle_left && lane > 0) {          
              lane--; 
            // case no vehicle right and it's possible to change
            } else if (!is_vehicle_right && lane != 2) {
              lane++; 
            } else {
              speed_diff -= .224;
            }
          } else {
            // case ego car not on center lane
            if (ref_vel < 49.5) {
              speed_diff += .224;
            }
          }
              
          // Create a list of evenly spaced (30m) waypoints(x,y)
          vector<double> ptsx;
          vector<double> ptsy;             

          // reference x, y, yaw 
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          // if previous size is almost empty, use the car as starting reference
          if(prev_size < 2) {
            // use two points that make the path tangent to the car
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);

            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          } else {
            // use previous path's end point as starting reference
            // redefine reference state as previous path end point
            ref_x = previous_path_x[prev_size-1];
            ref_y = previous_path_y[prev_size-1];

            double ref_x_prev = previous_path_x[prev_size-2];
            double ref_y_prev = previous_path_y[prev_size-2];
            ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

            // use two points that make the path tangent to the car
            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);

            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);  
          }

          // in frenet add evenly 30m spaced points ahead of the starting reference
          vector<double> next_wp0 = getXY(car_s+30,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s+60,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s+90,(2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_wp0[0]); 
          ptsx.push_back(next_wp1[0]); 
          ptsx.push_back(next_wp2[0]); 

          ptsy.push_back(next_wp0[1]); 
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);
                
          for (int i = 0; i < ptsx.size(); i++) {
            // shift car reference angle to 0 degree
            double shift_x = ptsx[i]-ref_x;
            double shift_y = ptsy[i]-ref_y;

            ptsx[i] = (shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
            ptsy[i] = (shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw)); 
          }
                
          // create a spline
          tk::spline s;
               
          // set (x,y) points to the spline
          s.set_points(ptsx, ptsy);

          // define the actual (x,y) points we will use for the planner
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // start with all of the previous path points from last time
          for(int i = 0; i < previous_path_x.size(); i++) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // calculate how to break up splines points so that we travel at our desired ref. velocity
          double target_x = 30.0;
          double target_y = s(target_x);
          double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

          double x_add_on = 0;
            
          // fill up the rest of our path planner after filling it with prev. points, here we output always 50 points
          for(int i=1; i <= 50-previous_path_x.size(); i++) {  
            //limit speed and acc
            ref_vel += speed_diff;
            if (ref_vel > 49.5) {
              ref_vel = 49.5;
            } else if (ref_vel < 0.224) {
              ref_vel = 0.224;
            }
            
            double N = (target_dist/(.02*ref_vel/2.24));
            double x_point = x_add_on+(target_x)/N;
            double y_point = s(x_point);

            x_add_on = x_point;

            double x_ref = x_point;
            double y_ref = y_point;

            // rotate back to normal after rotating it earlier
            x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
            y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));

            x_point += ref_x;
            y_point += ref_y;

            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);               
          }
          // END
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";
        
          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);   
          oldT = newT; 
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
