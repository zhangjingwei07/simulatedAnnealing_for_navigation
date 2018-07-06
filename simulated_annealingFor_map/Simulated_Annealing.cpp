/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Simulated_Annealing.h"

#include "vector"




//"Main function"
double Simulated_Annealing(vector< vector <double> >& route_time, vector<node_ID>& current_node_path){
    
    //set up a random instance
    
    TSP_Instance tsp_instance;
 
    
    tsp_instance.createNew(route_time, current_node_path.size(), current_node_path);
 
    
    //Initialize all related classes
    Optimizer optimizer;
    
    //Initialize three strategies of moving
    ChainReverseMove move1;
    SwapCityMove move2;
    RotateCityMove move3;
 
    
    //Add these strategies into optimize plan "moves" 
    optimizer.add_move(&move1);
    
    optimizer.add_move(&move2);
    
    optimizer.add_move(&move3);
    
    
    //Set up the parameter for annealing, including initial temperature, 
    //end temperature and delta, which is constant
    Geometric_Cooling_Action action(150, 1e-2, 0.95);
    
    //put the cooling schedule into optimizer
    optimizer.cooling_action = & action;
    
    //outer loop define the times to lower the temperature 
    optimizer.outer_loop = 100;
    
    //inner loop defines the 
    optimizer.inner_loop = 5000;
    

        
    //Start optimization
    double best = optimizer.Optimize(tsp_instance, current_node_path);   
    
    return best;
    
    
}


double TSP_Instance::calculate_tour_length(const std::vector<node_ID>& tour) const{
    
    assert(tour.size() == cities.size());
    
    double result = 0;

//********************************Do not consider first and last depot*********    
    for(size_t index = 0; index < tour.size() -1 ; index++){
        

        result = result + route_time[tour[index]][tour[index+1]];
    }
    //We need to consider the starting and ending depot
    result += globalVar.depot_distance[tour[0]];
    result += globalVar.depot_distance[tour[tour.size() - 1]];
    
    return double(result);
    
}


double time_difference_chain_reverse(vector< vector<double> > & route_time, std::vector<node_ID> & state, node_ID i, node_ID j){
   
    double time_difference = 0;
    
    if(i == 0){
        time_difference -= globalVar.depot_distance[state[i]];
        time_difference += globalVar.depot_distance[state[j]];
        
    }
    
    
    
    //calculate the conjunction point
    if(i != 0){
        time_difference -= route_time[state[i-1]][state[i]];
        time_difference += route_time[state[i-1]][state[j]];
        
    }
    
    //calculate the middle part
    for(node_ID d = i; d <= j - 1; d++){
        time_difference -= route_time[state[d]][state[d+1]];
        time_difference += route_time[state[d+1]][state[d]];
               
    }
    
    //if j is not the last index
    if(j <= state.size() -2 ){
        time_difference -= route_time[state[j]][state[j+1]];
        time_difference += route_time[state[i]][state[j+1]];
    }
       
    //else j is the end point
    if(j == state.size() - 1){
        time_difference -= globalVar.depot_distance[state[j]];
        time_difference += globalVar.depot_distance[state[i]];
    }
    
    
    return time_difference;
       
    
}



double time_difference_swap(vector< vector<double> > & route_time, std::vector<node_ID> & state, node_ID i, node_ID j){
    
    double time_difference = 0;
     
//    cout << "depot size: "<< globalVar.depot_distance.size() << endl;
//    cout << "state size: "<< state.size() << endl;
    
    //Take depot into consideration
    if(i == 0){
        time_difference -= globalVar.depot_distance[state[i]];
        time_difference += globalVar.depot_distance[state[j]];
        
    }
    
    //i is not starting point
    if(i != 0){
        time_difference -= route_time[state[i-1]][state[i]];
        time_difference += route_time[state[i-1]][state[j]];
        
    }
    
    
    //If they are adjacent, just swap them
    if(j == i+1){
        time_difference -= route_time[state[i]][state[j]];
        time_difference += route_time[state[j]][state[i]];
    }
    
   
    
    //They are not adjacent
    else{
        time_difference -= route_time[state[j-1]][state[j]];
        time_difference += route_time[state[j]][state[i+1]];
        
        time_difference -= route_time[state[i]][state[i+1]];
        time_difference += route_time[state[j-1]][state[i]];
        
    }
    
    if(j <= state.size()-2){
        time_difference -= route_time[state[j]][state[j+1]];
        time_difference += route_time[state[i]][state[j+1]];
    }
    
    //else j is the end point
    if(j == state.size() - 1){
        time_difference -= globalVar.depot_distance[state[j]];
        time_difference += globalVar.depot_distance[state[i]];
    }
    
    return time_difference;
 
}




double time_difference_rotate(vector< vector<double> > & route_time, std::vector<node_ID> & state, node_ID i, node_ID j, node_ID k){
    
    double time_difference = 0;
    
    //take the depot into consideration
    
    if(i == 0){
        
        time_difference -= globalVar.depot_distance[state[i]];
        time_difference += globalVar.depot_distance[state[j]];
    }
    
    if(i != 0){
        time_difference -= route_time[state[i-1]][state[i]];
        time_difference += route_time[state[i-1]][state[j]];
        
    }
            
    time_difference -= route_time[state[j-1]][state[j]];       
    
    if(k <= state.size()-2){
        time_difference -= route_time[state[k]][state[k+1]];
        time_difference += route_time[state[j-1]][state[k+1]];
    }
    
    //k is the end point, take depot into consideration
    else{
        time_difference -= globalVar.depot_distance[state[k]];
        time_difference += globalVar.depot_distance[state[j-1]];
        
    }
    
    time_difference += route_time[state[k]][state[i]];
 
    
    return time_difference;
}


//*****************The core function:: Optimizer********************************

double Optimizer::Optimize(const TSP_Instance & instance, std::vector<node_ID> & result) const{
    
    //First get the number of cities
    int n = static_cast<int>(instance.get_cities().size());
    
    assert(n > 0);
    
    assert(moves.size() > 0);
    
    
    //Set up the run time configuration
    Config config;
    
    
    
    //Set up the initial path
    config.state.resize(n);
    
    config.best_state.resize(n);
    
    
    //Put the node pattern into the state
    for(int index = 0; index < n;  index++){
        
        config.state[index] = instance.cities[index];
        config.best_state[index] = instance.cities[index];
        
    }
//    
//    cout <<"Initial solution" << endl;
//    cout << "state size:" << config.state.size() << endl;
//    
//    for(int i = 0; i < config.state.size(); i++)
//                cout << " "<< config.state[i];
//    
//    cout << endl;
    
    
    
    
    //Calculate the initial energy, which is the tour length
    //Remember to add the start and end depot
    config.energy = instance.calculate_tour_length(config.state);

    
    
    
    
    //cout << "Current Energy: " << config.energy << endl;
   
     
    
    
    config.best_energy = config.energy;
    
    //Set up the initial temperature
    config.temp = cooling_action->initial_temp();
    
    
    
    std::mt19937 g({std::random_device{}()});
    
    //set up an initial distribution over the possible moves, range from 0 to moves index
    std::uniform_int_distribution<int> moveDist(0, static_cast<int> (moves.size() -1));
    
    //A uniform distribution for the acceptance probability, range from 0 to 1  
    std::uniform_real_distribution<double> uniformDist(0.0, 1.0);
    
    
    //set up the mover service, n means the number of nodes
    Optimizer::Move_Service* service = new Optimizer::Move_Service(n);   
    
    
    for(size_t i = 0; i < moves.size(); i++){
        
        moves[i]->set_Move_Service(service);
        
    }
    
    //The current proposal
    std::vector<node_ID> proposal;
    
    

    
    //start the optimization
    for(config.out_loop = 0; config.out_loop < outer_loop; config.out_loop++){
        
        //Determine the next temperature, if smaller than end........
        config.temp = cooling_action-> next_temp(config);
        
        
        //Simulate the process of Markov Chain
        for(config.inner_loop = 0; config.inner_loop < inner_loop; config.inner_loop++){
            
            //proposal means the current result, which is not optimal
            proposal = config.state;
            
            
            
            //propose a new neighbor according to some moves
            //choose a move
            int m = moveDist(g);
          
            
            //return the path difference
            const double delta = moves[m]->propose(proposal);
            
      //      cout << "Energy1: " << config.energy << endl;
           
            
            //Get the energy of the new proposal
            const double energy = config.energy + delta; 

            double result1 = 0;
            
     //       cout << "Energy2: " << energy << endl;
            
            
          //  result1 = instance.calculate_tour_length(config.state);
  
            
         //   cout << "using formola new energy: " << result1 << endl;
//            
//        
//        
//        
//            cout << "New_energy: " << energy << endl;
//            
//            cout << "Delta: " << delta << endl;
          
            
            /*******************************/
            //assert(energy > 0);
             /***************************/
            
            //It means that the new path takes less time
            if (delta <= 0){
                
                //Accept the move
                config.state = proposal;
                
                config.energy = energy;
            
                if(energy < config.best_energy){
                    
                    //This energy is better
                    config.best_energy = energy;
                    config.best_state = proposal;
                       
                    //Used to debug
//                   if (config.best_state == pattern) {
//                       bool x = 0;
//                       assert(x == 1);
                                          
                   //}
                    
                }
            
            }
             else  {
                
                
                //Accept the proposal with certain probability
                
                double u = uniformDist(g);
                
//                cout << "random generator: " << u << endl;
//                cout << "The current threshold: " << exp(-1/config.temp * delta) << endl;
//                cout << "current temp: "<< config.temp << endl;
                if(u <= std::exp(-1/config.temp * delta)){
                    
                    config.state = proposal;
                    config.energy = energy;
                    
                    
                }
                
             }
            
            
            
            
            
//            cout << "solution : " << endl;
//            
//            for(int i = 0; i < config.state.size(); i++){
//                cout << "  "<< config.state[i];
//            }
//            
//            cout << endl;
          
        
        }
            
        
        
        auto current = std::chrono::high_resolution_clock::now();    
            
        auto alpha = std::chrono::duration_cast<chrono::duration<double>>(current - globalVar.STARTING_TIME);
 
        
        
//        cout << "time now: " << alpha.count() << endl;
        
        if(alpha.count() >= 0.95 * TIME_LIMIT)
          
            return config.best_energy;
         
    }
    
    DELETE_PTR(service);
    
    for(size_t i = 0; i < moves.size(); i++){
        
        moves[i]-> set_Move_Service(0);
    }
    
    result = config.best_state;
    
    
    //Do the final notification
    config.finish = true;
    config.state = config.best_state;
    config.best_energy = config.best_energy;
  
    
  //  cout << "BEST ENERGY:" << config.best_energy << endl;
//    cout << "Final tour length" << instance.calculate_tour_length(config.best_state) << endl;
    
//    for(node_ID index = 0; index < config.best_state.size(); index++){
//        cout << config.best_state[index] << " ";
//        
//    }
   // cout << endl;
    return config.best_energy;
}
