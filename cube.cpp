
/*

*********************************
Project LSH
@Author Nikolas Gialitsis
AM:1115201400027
*********************************

*/

#include <iostream>
#include <string>
#include <cstdlib> // for exit
#include <fstream>
#include <stdlib.h>
#include <vector> 
#include <random>
#include "lsh_classes.h"
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <numeric>//std::inner_product
#include <unordered_map>
#include <math.h>
#include <complex> // std::complex, std::norm
#include <ctime>
#include <bitset>
#include <chrono>
using namespace std;



int modulo(int x,int y){
	return ( (x%y)+y)%y;

}

string euclidean_f(std::vector<int> h){
	string s;
	s.clear();
	
	for(int i = 0 ; i < h.size() ; i++){

		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::uniform_real_distribution<data_type> uniform_distribution(0.0,1.0);
		
		data_type fi = uniform_distribution(generator);
 		std::uniform_int_distribution<> dis(0,1);
 		long long random_num = floor(10*fi + h.at(i));

		if((random_num % 2) == 0)
			s.push_back('0');
		else
			s.push_back('1'); 
	}
	return s;

}

void EuclideanCube(int n,int K,int d,int W,string query_path,
	std::vector<Point>& dataset,int M,int probes,string& output_path){
	
	int vertices_num = pow(2,K);
	std::vector<Point> Cube[vertices_num];
	



	data_type e = 0.05;
	std::vector<data_type> v[K];
	space += (K*sizeof(data_type));
	data_type t[K];
	space += (K*sizeof(data_type));

  	std::random_device normal_generator;
	normal_distribution<double> normal_distribution(0.0,1.0);	

	std::random_device generator;
	uniform_real_distribution<double> uni_distribution(0.0,W-1);

	if(probes > vertices_num)probes=vertices_num;

	for(int k = 0 ; k < K ; k++){
		v[k].clear();
		for(int j = 0 ; j  < d ; j++){
			double item = normal_distribution(normal_generator);
			v[k].push_back(item);
		}

		t[k] = uni_distribution(generator);
	}
	
	for(int i = 0 ; i < n ; i++){

		std::cout<<"(euclidean)point "<<i<<std::endl;	
		int bucket = -1;

		std::vector<int> h;
		h.clear();

		for(int k = 0 ; k < K ; k++){
			int h_val = euclidean_h(dataset.at(i),W,d,v[k],t[k]);
			h.push_back(h_val);
		}
		assert(h.size() == K);
		
		
		string gs;
		gs.assign(euclidean_f(h));
		assert(gs.size() == K);
		int verticle = BinaryToDecimal(gs);
		assert(verticle >= 0);

		Cube[verticle].push_back(dataset.at(i));
	
	}

	std::cout<<"----- Printing Cube -------"<<std::endl<<std::endl;

	for(int verticle = 0 ; verticle < vertices_num; verticle++){
		int vertices = 0;
		std::cout <<"\tvertice: "<< verticle << " : ";
		for(vector<Point>::iterator it = Cube[verticle].begin() ; 
			it != Cube[verticle].end() ;it++ ){
				vertices++;
		}
		std::cout<<vertices<<std::endl;
	}



	data_type R = 350;//default

	vector<Point> queries; 
	ifstream input_queries;
	int query_points = 0;
	input_queries.open(query_path.c_str());

	bool input_radius = false;
	if (input_queries.is_open()){
		string line;
		while ( getline (input_queries,line)){

			if(line.size() == 1)continue;

			if(query_points == 0){
				string radius_v;
				for(int j = 0 ; j < line.size() ; j++){					
					if(line[j] != '\r'){
						radius_v.push_back(line[j]);
					}
					if((line[j] == ' ')||(line[j] == '\t')){
						string eq;
						eq.assign("Radius: ");
						if(radius_v == eq){
							radius_v.clear();
							input_radius = true;						
						}
					}
					if(j == (line.size()-1)){
						if(input_radius){							
							R = atof(radius_v.c_str());
							break;
						}
					}
				}
			}
			if(input_radius){
				input_radius = false;
				continue;
			}
			Point* ptr = new Point(line);
			queries.push_back(*ptr);
			query_points++;
		}

		
	}
	else cout << "Unable to open file"<<endl;
	

	cout<<"R = |"<<R<<"|"<<endl;


	ofstream output;
	output.open(output_path.c_str());
	if (!output.is_open()){
		std::cout<<"Output path not correct"<<endl;
		exit(-2);
	}



	int found_items = 0;
	float average_time = 0;
	data_type approx = 0;
	data_type max_ratio = -1;

	
	for(int query = 0 ; query < query_points; query++){
		
		output<<std::endl<<"--------------------"<<std::endl;
		output<< "********************"<<std::endl;
		output<<"--------------------"<<std::endl;
		output<<"(euclidean)Query:"<<query<<std::endl;
		std::cout<<"(euclidean)Query:"<<query<<std::endl;
		const clock_t begin_time = clock();

		Point p = queries.at(query); 	
		

		int bucket = -1;

		std::vector<int> h;
		h.clear();

		for(int k = 0 ; k < K ; k++){
			int h_val = euclidean_h(p,W,d,v[k],t[k]);
			h.push_back(h_val);
		}
		assert(h.size() == K);
		space += (K*sizeof(int));
		
		string gs;
		gs.assign(euclidean_f(h));
		assert(gs.size() == K);
		space += (K*sizeof(char));


		int verticle = BinaryToDecimal(gs);
		assert(verticle >= 0);
		
		int hamming_lvl = 0;

		bool visited [vertices_num];
		for(int i = 0 ; i < vertices_num ;i++)
			visited[i] = false;
		space += (vertices_num*sizeof(bool));
		data_type min_dist = 10*R;


		int num_items = 0;
		int nearest_id = -1;
		for(int counter = 0 ; counter < probes ; counter++){

			if(hamming_lvl == 0)
				hamming_lvl++;
			else{

				bool found_verticle =false;
				while(found_verticle != true){

					for(int neighbour = 0 ; neighbour < vertices_num ; neighbour++){
						string large_bin =  std::bitset<100>(neighbour).to_string();
						string  ns;
						ns.clear();

						for(int s = 100-K;s < 100 ; s++){
							ns.push_back(large_bin.at(s));
						}
						space += (K*sizeof(char));

						if(hammingDistance(gs,ns) == hamming_lvl){
							if(visited[neighbour] == false){
								verticle = neighbour;
								visited[neighbour] = true;
								found_verticle = true;
								break;
							}
						}
					}
					if(found_verticle == false)
						hamming_lvl++;
				}

			}
			for(vector<Point>::iterator it = Cube[verticle].begin() ; 
				it != Cube[verticle].end() ;it++ ){
				num_items++;

				data_type dist = EuclideanDistance(*it,p);
				assert(dist >= 0);
				if(dist <= (1+e)*R){			
					num_items++;
					if(num_items > M)break;
					output<<"Item"<<it->get_id()<<std::endl;
					if(dist < min_dist){
						min_dist = dist;
						nearest_id = it->get_id();
					}

				}
				
			}		
		}

		




		float tCube =  (float( clock () - begin_time ) /  CLOCKS_PER_SEC);
		

		

		average_time += tCube;

		data_type true_min_dist=10*R;


		const clock_t begin_time_2 = clock();
		for(int i = 0 ; i < n; i++){
			data_type dist = EuclideanDistance(p,dataset.at(i));
			if(dist < true_min_dist)true_min_dist = dist;
		}

		output<<"Nearest neighbour: Item  ";
		if(nearest_id != -1){
			output<<nearest_id<<std::endl;
			output<<"distanceCube :" << min_dist << std::endl;
		}
		else{
			output<<"not found"<<std::endl;
		}

		output<<"distanceTRUE :" << true_min_dist << std::endl;
		output<<"tCube: "<<tCube<<" sec"<<endl;
		output<<"tTrue:"<<(float( clock () - begin_time_2 ) /  CLOCKS_PER_SEC)<<" sec"<<endl;
		


		assert(min_dist >= 0);
		assert(true_min_dist >= 0);
		assert(min_dist >= true_min_dist);
		if(nearest_id != -1){
			data_type p_error = min_dist/true_min_dist;
			if(p_error >= max_ratio)max_ratio = p_error;
			approx += p_error;
			found_items++;
		}


		
	}


	output<<"*********************"<<std::endl;
	output<<"Average Query Time : "<<average_time/((float)query_points) << " seconds"<<std::endl;
	output<<"Approx "<<approx/(float)found_items<<std::endl;
	output<<"Max Ratio :"<<max_ratio<<std::endl;
	output<<"Space :"<<space<<" bits"<<std::endl;
	output<<"*********************"<<std::endl;


	input_queries.close();

}


unsigned b2d(unsigned num)
{
    unsigned res = 0;

    for(int i = 0; num > 0; ++i)
    {
        if((num % 10) == 1)
            res += (1 << i);

        num /= 10;
    }

    return res;
}


int BinaryToDecimal(string& h){
	return b2d(atoi(h.c_str()));

} 


void euclidean_g(std::vector<int>& gv,std::vector<int>& h,std::vector<int>& shuffle,int K){
	
	gv.clear();
	for(int i = 0 ; i  < K ; i++){
		int ident = shuffle[i];
		assert(ident >= 0);
		//gv.push_back(h.at(i));
		gv.push_back(h.at(shuffle[i]));
	}

} 


double EuclideanDistance(Point& p1,Point& p2){


	string s1;
	s1.clear();
	s1.assign(p1.get_all_fields());
	boost::tokenizer<> tok(s1);

	std::vector<data_type> v1;
	v1.clear();
	for (boost::tokenizer<>::iterator it=tok.begin(); it!=tok.end();it++) {
		string s;
		s.assign(*it);
		data_type item = atof(s.c_str());
		v1.push_back(item);
	}

	string s2;
	s2.clear();
	s2.assign(p2.get_all_fields());
	boost::tokenizer<> tok2(s2);
	std::vector<data_type> v2;
	v2.clear();
	for (boost::tokenizer<>::iterator it=tok2.begin(); it!=tok2.end();it++) {;
		string s;
		s.assign(*it);
		data_type item = atof(s.c_str());
		v2.push_back(item);
	}


	data_type dist = 0 ;

	assert(v1.size() == v2.size());

	for(int i = 0 ; i < v1.size(); i++){
		data_type coord_dist = v1.at(i)-v2.at(i);
		dist = dist + coord_dist*coord_dist;// pow(coord_dist,2);
	}
	if(dist != 0)
		dist = sqrt(dist);
	assert(dist >= 0);

	
	return dist;
}


int euclidean_h(Point& point,int w,int d,std::vector<data_type>& v,data_type t){



	string s;
	s.assign(point.get_all_fields());
	boost::tokenizer<> tok(s);

	std::vector<data_type> p;
	p.clear();
	for (boost::tokenizer<>::iterator it=tok.begin(); it!=tok.end();it++) {
		data_type item = atof(it->c_str());
		p.push_back(item);
	}


	double init_value = 0;
	double pv = (data_type)inner_product(p.begin(),p.begin()+p.size(),v.begin(),init_value);
	


	assert(p.size() == d);
	assert(v.size() == p.size());
	assert(w > 0);
	double result = (pv+t)/(double(w));
	return (int)result;
}


double CosineDistance(Point& p1,Point& p2){
	

	string s1 = p1.get_all_fields();
	boost::tokenizer<> tok(s1);

	std::vector<data_type> v1;
	v1.clear();
	for (boost::tokenizer<>::iterator it=tok.begin(); it!=tok.end();it++) {
		data_type item = atof(it->c_str());
		v1.push_back(item);
	}

	string s2 = p2.get_all_fields();
	boost::tokenizer<> tok2(s2);
	std::vector<data_type> v2;
	v2.clear();
	for (boost::tokenizer<>::iterator it=tok2.begin(); it!=tok2.end();it++) {;
		data_type item = atof(it->c_str());
		v2.push_back(item);
	}

	data_type init_value = 0;
	assert(v1.size() == v2.size());
	data_type xy = (data_type)inner_product(v1.begin(),v1.begin()+v1.size(),v2.begin(),init_value);


	double sum_v1 = 0.0;
	double sum_v2 = 0.0;

	for(int i = 0 ; i < v1.size() ; i++){
		sum_v1 = sum_v1 + v1.at(i)*v1.at(i);
		sum_v2 = sum_v2 + v2.at(i)*v2.at(i);
	}
	data_type similarity = (xy)/(sqrt(sum_v1) * sqrt(sum_v2));
	data_type dist = 1 - similarity;
	return dist;
}

int hammingDistance(string& s1,string& s2){
	assert(s1.size() == s2.size());
	int off_bits = 0;
	for(int i = 0 ; i < s1.size() ; i++){
		if(s1.at(i) != s2.at(i))
			off_bits++;
	}

	return off_bits;

}


void CosineCube(int n,std::vector<Point>& dataset,int K,int d,
	string& query_path,int M,int probes,string& output_path){

	double e = 0.05;
	vector<data_type> r[K];
	space += (K*sizeof(data_type));

	std::random_device normal_generator;
	normal_distribution<double> normal_distribution(0.0,1.0);

	for(int k = 0 ; k < K; k++){ //plithos h	
		for(int y = 0 ; y < d ; y++){
			double item = normal_distribution(normal_generator);
			r[k].push_back(item);
		}
	}
	
		
	int vertices_num = pow(2,K);
	std::vector<Point> Cube[vertices_num];


	for(int i = 0 ; i < n ; i++){
		std::cout<<"(cosine) point "<<i<<std::endl;
	
		string gs;			
		for(int k = 0 ; k < K  ; k++){
			gs.push_back(cosine_h(dataset.at(i),r[k],K));
		}	

		int verticle = BinaryToDecimal(gs);
		assert(verticle >= 0);
		Cube[verticle].push_back(dataset.at(i));
			
		
	}
	
	std::cout<<"----- Printing Cube -------"<<std::endl<<std::endl;

	for(int verticle = 0 ; verticle < vertices_num; verticle++){
		int vertices = 0;

		std::cout <<"\tvertice: "<< verticle << " : ";
		for(vector<Point>::iterator it = Cube[verticle].begin() ; 
			it != Cube[verticle].end() ;it++ ){
				vertices++;
		}
		std::cout<<vertices<<std::endl;
	}

	//Parsing output file - store all items in a vector



	data_type R = 1.0;//default

	vector<Point> queries; 
	ifstream input_queries;
	int query_points = 0;
	input_queries.open(query_path.c_str());

	bool input_radius = false;
	if (input_queries.is_open()){
		string line;
		while ( getline (input_queries,line)){

			if(line.size() == 1)continue;

			if(query_points == 0){
				string radius_v;
				for(int j = 0 ; j < line.size() ; j++){					
					if(line[j] != '\r'){
						radius_v.push_back(line[j]);
					}
					if((line[j] == ' ')||(line[j] == '\t')){
						string eq;
						eq.assign("Radius: ");
						if(radius_v == eq){
							radius_v.clear();
							input_radius = true;						
						}
					}
					if(j == (line.size()-1)){
						if(input_radius){							
							R = atof(radius_v.c_str());
							break;
						}
					}
				}
			}
			if(input_radius){
				input_radius = false;
				continue;
			}
			Point* ptr = new Point(line);
			queries.push_back(*ptr);
			query_points++;
		}

		
	}
	else cout << "Unable to open file"<<endl;
	

	cout<<"R = |"<<R<<"|"<<endl;

	ofstream output;
	output.open(output_path.c_str());
	if (!output.is_open()){
		cout<<"Output path not correct"<<endl;
		exit(-1);
	}



	double  average_time = 0;
	double approx = 0;
	int found_items = 0;
	if(probes > vertices_num)probes=vertices_num;
	data_type max_ratio = -1;


	
	for(int query = 0 ; query < query_points; query++){

		output<<std::endl<<"--------------------"<<std::endl;
		output<<           "********************"<<std::endl;
		output<<"--------------------"<<std::endl;

		output<<"(cosine)Query:"<<query<<std::endl;
		std::cout<<"(cosine)Query:"<<query<<std::endl;

		data_type min_dist = 10*R;
		const clock_t begin_time = clock();
		
		Point p = queries.at(query); 
		int items_found = 0;
	
		string gs;	
		for(int k = 0 ; k < K  ; k++){
			gs.push_back(cosine_h(p,r[k],K));
		}

		int verticle = BinaryToDecimal(gs);
		assert(verticle >= 0);
		int num_items = 0;
		int hamming_lvl = 0;
		int nearest_id = -1;
		bool visited [vertices_num];
		space += (vertices_num*sizeof(bool));
		for(int i = 0 ; i < vertices_num ;i++)visited[i] = false;


		
		output<<"R-near neighbours:"<<std::endl;
		for(int counter = 0 ; counter < probes ; counter++){

			if(hamming_lvl == 0)
				hamming_lvl++;
			else{

				bool found_verticle =false;
				while(found_verticle != true){
					for(int neighbour = 0 ; neighbour < vertices_num ; neighbour++){
						string large_bin =  std::bitset<100>(neighbour).to_string();
						string  ns;
						ns.clear();

						for(int s = 100-K;s < 100 ; s++){
							ns.push_back(large_bin.at(s));
						}
						space += (K*sizeof(char));
						if(hammingDistance(gs,ns) == hamming_lvl){
							if(visited[neighbour] == false){
								verticle = neighbour;
								visited[neighbour] = true;
								found_verticle = true;
								break;
							}
						}
					}
					if(found_verticle == false)
						hamming_lvl++;
				}

			}



			for(vector<Point>::iterator it = Cube[verticle].begin() ; 
				it != Cube[verticle].end() ;it++ ){
				

				data_type dist = CosineDistance(*it,p);
				assert(dist >= 0);
				if(dist <= (1+e)*R){			
					num_items++;
					if(num_items > M)break;

					output<<"Item"<<it->get_id()<<std::endl;
					if(dist <= min_dist){
						min_dist = dist;
						nearest_id = it->get_id();
					}

				}
				
			}		
		}

		
		

		float tCube =  (float( clock () - begin_time ) /  CLOCKS_PER_SEC);
		
		average_time += tCube;

		data_type true_min_dist=10*R;

	
		const clock_t begin_time_2 = clock();		
		for(int i = 0 ; i < n; i++){
			data_type dist = CosineDistance(p,dataset.at(i));
			if(dist < true_min_dist)true_min_dist = dist;
		}


		output<<"Nearest neighbour: Item  ";
		if(nearest_id != -1){
			output<<nearest_id<<std::endl;
			output<<"distanceCube :" << min_dist << std::endl;
		}
		else{
			output<<"not found"<<std::endl;
		}

		output<<"distanceTRUE :" << true_min_dist << std::endl;
		output<<"tCube: "<<tCube<<" sec"<<endl;
		output<<"tTrue:"<<(float( clock () - begin_time_2 ) /  CLOCKS_PER_SEC)<<" sec"<<endl;
		
		
		


		assert(min_dist >= true_min_dist);
		if(nearest_id != -1){
			data_type p_error = min_dist/true_min_dist;
			if(p_error >= max_ratio)max_ratio = p_error;
			approx += p_error;
			found_items++;
		}


		
	}

	output<<"*********************"<<std::endl;
	output<<"Average Query Time : "<<average_time/((float)query_points) << " seconds"<<std::endl;
	output<<"Approx "<<approx/(float)found_items<<std::endl;
	output<<"Max Ratio :"<<max_ratio<<std::endl;
	output<<"Space :"<<space<<" bits"<<std::endl;
	output<<"*********************"<<std::endl;

	input_queries.close();


}



char cosine_h(Point& point,vector<data_type>& r,int K){
	
	string s;
	s.assign(point.get_all_fields());
	boost::tokenizer<> tok(s);

	std::vector<data_type> p;
	p.clear();
	for (boost::tokenizer<>::iterator it=tok.begin(); it!=tok.end();it++) {
		string s;
		s.assign(*it);
		data_type item = atof(s.c_str());
		p.push_back(item);
	}

	if((data_type)inner_product(p.begin(),p.begin()+p.size(),r.begin(),0) >= 0){
		return '1';
	}
	else{
		return '0';
	}	
}



int main(int argc,char *argv[]){

	if ((argc < 3)and(argc > 1)){
		cout << "Too few input arguments to function.(argc = "<<argc<<", expected > 3)"<<endl;
		exit(0);
	}


	string input_path;
	string query_path;
	string output_path;
	string metric;
	int K = 3;
	int d = 0;
	int n;
	int M = 10;

	const int W = 4;
	int probes = 2;
	//Initialization
	if(argc == 1){
		input_path.assign("./input_small");
		query_path.assign("./query_small");
	}

	for(int i = 1 ; i < argc ; i++){
		string s(argv[i]);
		if((s == "-d")and(argc > i+1)){
			string p = argv[i+1];
			input_path.assign(p);
		}
		else if((s == "-q")and(argc > i+1)){
			string p = argv[i+1];
			query_path.assign(p);
		}
		else if((s == "-k")and(argc > i+1)){
			K = atoi(argv[i+1]);
		}
		else if((s == "-o")and(argc > i+1)){
			string p = argv[i+1];
			output_path.assign(p);
		}
		else if((s == "-M")and(argc > i+1)){
			M = atoi(argv[i+1]);
		}
		else if((s == "-probes")and(argc > i+1)){
			probes = atoi(argv[i+1]);
		}		
		else;

	}

	vector<Point> dataset; 
	ifstream input;
	n  = -1;
	input.open(input_path.c_str());
	if (input.is_open()){
		string line;
		while ( getline (input,line)){
			if(n == -1){
				bool delimeter = false;
				if(line.size() == 0){
					metric.assign("euclidean");
				}
				for(int j = 0 ; j < line.size() ; j++){					

					if(line[j] != '\r'){
						metric.push_back(line[j]);
						if((line[j] == ' ')||(line[j] == '\t')){							
							delimeter = true;
							break;
						}
					}
				}
				if((delimeter == true)||(line.size() == 1))
					metric.assign("euclidean");
				n++;
				continue;
			}

			Point* ptr = new Point(line);
			dataset.push_back(*ptr);
			if( n == 0){
				string s;
				s.assign(ptr->get_all_fields());
				boost::tokenizer<> tok(s);
				for (boost::tokenizer<>::iterator it=tok.begin(); it!=tok.end();it++) {					
					data_type item = atof(it->c_str());
					d++;
				}
			}
			n++;
		}

		input.close();
	}
	else cout << "Unable to open file"<<endl; 


	

	if(metric == "euclidean")
		EuclideanCube(n,K,d,W,query_path,dataset,M,probes,output_path);
	else if(metric == "cosine")
		CosineCube(n,dataset,K,d,query_path,M,probes,output_path);
	else
		EuclideanCube(n,K,d,W,query_path,dataset,M,probes,output_path);
	return 0;



}