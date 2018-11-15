
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

using namespace std;

int long_modulo(int x,long long y){
	return ( (x%y)+y)%y;

}

int modulo(int x,int y){
	return ( (x%y)+y)%y;
}


void EuclideanLSH(int n,int L,int K,int d,int W,string query_path,std::vector<Point>& dataset,string& output_path ){
	
	int TableSize = n/4;
	data_type e = 0.05;
	vector<Point> Hash_Tables[L][TableSize];
	std::vector<data_type> v[L][K];
	data_type t[L][K];
	space += (L*K*sizeof(data_type));

  	std::random_device normal_generator;
	normal_distribution<double> normal_distribution(0.0,1.0);	

	std::random_device generator;
	uniform_real_distribution<double> uni_distribution(0.0,W-1);

	for(int l = 0 ; l < L ; l++){
		
		for(int k = 0 ; k < K ; k++){
			v[l][k].clear();
			for(int j = 0 ; j  < d ; j++){
				double item = normal_distribution(normal_generator);
				v[l][k].push_back(item);
			}

			t[l][k] = uni_distribution(generator);


		}
		
	}

	std::random_device rd;
	uniform_int_distribution<int> uniform_distribution(-20,20);
	std::vector<int> r[L];	
	for(int i = 0 ; i < L ; i++){
		r[i].clear();		
		for(int j = 0 ; j < K; j++){
			int item = uniform_distribution(rd);
			r[i].push_back(item);
		}
	}
	space += (L*K*sizeof(int));

	std::vector<int> shuffle_idents[L];
	space += (L*sizeof(int));
	std::random_device shuffle_generator;
	uniform_int_distribution<int> shuffle_distribution(0,K-1);

	bool diff=false;
	while(diff != true){
		for(int i = 0 ; i < L ; i++){
			for(int j = 0 ; j < K;j++){
				int item = shuffle_distribution(shuffle_generator);
				shuffle_idents[i].push_back(item);
			}
			diff = true;
			for(int k = 0 ; k < i ; k++){
				if(shuffle_idents[i] == shuffle_idents[k]){
					diff = false;
					break;
				}
			}
			if(diff == false){
				break;
			}
		}

	}

	for(int i = 0 ; i < n ; i++){

		std::cout<<"(euclidean)point "<<i<<std::endl;
		
		int bucket = -1;

		std::vector<int> gv[L];
		
		

		for(int l = 0 ; l < L ; l++){
			
			std::vector<int> h;
			h.clear();

			for(int k = 0 ; k < K ; k++){
				int h_val = euclidean_h(dataset.at(i),W,d,v[l][k],t[l][k]);
				h.push_back(h_val);
			}
			assert(h.size() == K);


			gv[l].clear();
			euclidean_g(gv[l],h,shuffle_idents[l],K);

		}

		for(int k = 0 ; k < L ; k++){
			int bucket = euclidean_phi(r[k],gv[k],K,TableSize);
			vector<Point>* myv = &(Hash_Tables[k][bucket]);			
			dataset.at(i).g = gv[k];
			space += (gv[k].size()*sizeof(char));
			myv->push_back(dataset.at(i));
		

		}
	}



	std::cout<<"----- Printing Hash -------"<<std::endl<<std::endl;

	for(int h = 0 ; h < L ; h++){
		std::cout <<"Hash Table "<< h << " : "<< std::endl;
		
		for(int b = 0 ; b < TableSize; b++){
			int bucket_num = 0;
			std::cout <<"\tBucket: "<< b << " : ";


			for(vector<Point>::iterator it = Hash_Tables[h][b].begin() ; 
				it != Hash_Tables[h][b].end() ;it++ ){
				bucket_num++;
			}
			std::cout<<bucket_num<<std::endl;

		}

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
	output.open(output_path.c_str(),std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open()){
		cout<<"Output path not correct"<<endl;
		exit(-1);
	}


	int found_items = 0;
	float average_time = 0;
	data_type approx = 0;

	data_type max_ratio = -1;

	
	for(int query = 0 ; query < query_points; query++){
		
		output<<std::endl<<"--------------------"<<std::endl;
		output<<std::endl<<"--------------------"<<std::endl;
		output<<"(euclidean)Query:"<<query<<std::endl;
		std::cout<<"(euclidean)Query:"<<query<<std::endl;
		const clock_t begin_time = clock();

		Point p = queries.at(query); 
	
		
		std::vector<int> gv[L];	

		
		output<<"R-near neighbours:"<<std::endl;
		for(int l = 0 ; l < L ; l++){

			std::vector<int> h;
			h.clear();
			for(int k = 0 ; k < K ; k++){
				int h_val = euclidean_h(p,W,d,v[l][k],t[l][k]);
				h.push_back(h_val);
			}
			assert(h.size() == K);



			gv[l].clear();
			euclidean_g(gv[l],h,shuffle_idents[l],K);							
		}

		data_type true_min_dist = 10*R;
		int bucket = -1;
		int items_found = 0;
		data_type min_dist = 10*R;
		int nearest_id = -1;

		for(int k = 0 ; k < L ; k++){
			int bucket = euclidean_phi(r[k],gv[k],K,TableSize);
			assert(bucket >= 0);				
			for(vector<Point>::iterator it = (Hash_Tables[k][bucket]).begin() ; 
				it != (Hash_Tables[k][bucket]).end() ;it++ ){

				if(it->g != gv[k]){
					continue;
				}

				data_type dist = EuclideanDistance(*it,p);
				assert(dist >= 0);
				if(dist <= (1+e)*R){			

					items_found++;
					output<<"Item"<<it->get_id()<<std::endl;

					if(dist < min_dist){
						min_dist = dist;
						nearest_id = it->get_id();
					}

					if(items_found > 3*L)break;

				}
				
			}	


			if(nearest_id == -1){
				if(bucket == TableSize-1)
					bucket = TableSize-2;
				else if(bucket == 0)
					bucket = 1;
				else{
					std::random_device coin_toss_generator;
					uniform_int_distribution<int> coin_distribution(0,1);
					int item = coin_distribution(coin_toss_generator);
					if(item == 0)
						bucket = bucket - 1;					
					else
						bucket = bucket + 1;
				}
				for(vector<Point>::iterator it = (Hash_Tables[k][bucket]).begin() ; 
					it != (Hash_Tables[k][bucket]).end() ;it++ ){

					if(it->g != gv[k]){
						continue;
					}

					data_type dist = EuclideanDistance(*it,p);
					assert(dist >= 0);
					if(dist <= (1+e)*R){			

						items_found++;
						output<<"Item"<<it->get_id()<<std::endl;

						if(dist < min_dist){
							min_dist = dist;
							nearest_id = it->get_id();
						}

						if(items_found > 3*L)break;
					}					
				}	
			}	
		}


		
		float tLSH =  (float( clock () - begin_time ) /  CLOCKS_PER_SEC);
		average_time += tLSH;



		const clock_t begin_time_2 = clock();
		for(int i = 0 ; i < n; i++){
			data_type dist = EuclideanDistance(p,dataset.at(i));
			if(dist < true_min_dist)true_min_dist = dist;
		}

		output<<"Nearest neighbour: Item  ";
		if(nearest_id != -1){
			output<<nearest_id<<std::endl;
			output<<"distanceLSH :" << min_dist << std::endl;
		}
		else{
			output<<"not found"<<std::endl;
		}

		output<<"distanceTRUE :" << true_min_dist << std::endl;
		output<<"tLSH: "<<tLSH<<" sec"<<endl;
		output<<"tTrue:"<<(float( clock () - begin_time_2 ) /  CLOCKS_PER_SEC)<<" sec"<<endl;
		

		assert(min_dist >= true_min_dist);
		
		if(nearest_id != -1){
			data_type p_error = min_dist/true_min_dist;
			if(p_error >= max_ratio)max_ratio=p_error;
			approx += p_error;
			found_items++;
		}
	}


	output<<"*********************"<<std::endl;
	output<<"Average Query Time : "<<average_time/((float)query_points) << " seconds"<<std::endl;
	output<<"Approx "<<approx/found_items<<std::endl;
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


int cosine_g(string& h,std::vector<int>& shuffle,int K){

	string result;
	result.clear();

	for(int i = 0 ; i  < K ; i++){
		int ident = shuffle[i];
		assert(ident >= 0);
		result.push_back(h.at(shuffle[i]));
	}

	return b2d(atoi(result.c_str()));

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



int euclidean_phi(std::vector<int>& r,std::vector<int>& g,int K,int TableSize){

	assert(g.size() == r.size());
	assert(g.size() > 0);
	long long rh = inner_product(g.begin(),g.end()+g.size(),r.begin(),0);


	long long M = pow(2,52)-5;
	assert(TableSize > 0);
	int temp_mod = long_modulo(rh,M);
	int bucket = modulo(temp_mod,TableSize);
	//int bucket = modulo(modulo(rh,M),TableSize);
	assert(bucket >= 0);
	return bucket;


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



void CosineLSH(int n,std::vector<Point>& dataset,int K,int L,int d,string& query_path,string& output_path){

	int TableSize = pow(2,K);
	double e = 0.05;
	std::cout<<"TableSize = "<<TableSize<<std::endl;
	assert(TableSize > 0);
	assert(L > 0);
	vector<Point> Hash_Tables[L][TableSize];

	vector<data_type> r[L][K];

	std::random_device normal_generator;
	normal_distribution<double> normal_distribution(0.0,1.0);


	for(int l = 0 ; l < L ; l++){
		for(int k = 0 ; k < K; k++){ //plithos h	
			for(int y = 0 ; y < d ; y++){
				double item = normal_distribution(normal_generator);
				r[l][k].push_back(item);
			}
		}
	}
		


	
	std::vector<int> shuffle_idents[L];
	std::random_device shuffle_generator;
	uniform_int_distribution<int> shuffle_distribution(0,K-1);

	bool diff=false;
	while(diff != true){
		for(int i = 0 ; i < L ; i++){
			for(int j = 0 ; j < K;j++){
				int item = shuffle_distribution(shuffle_generator);
				shuffle_idents[i].push_back(item);
			}
			diff = true;
			for(int k = 0 ; k < i ; k++){
				if(shuffle_idents[i] == shuffle_idents[k]){
					diff = false;
					break;
				}
			}
			if(diff == false){
				break;
			}
		}

	}


	for(int i = 0 ; i < n ; i++){
		std::cout<<"(cosine) point "<<i<<std::endl;

		for(int l = 0 ; l < L ; l++){
			string gs;			
			for(int k = 0 ; k < K  ; k++){
				gs.push_back(cosine_h(dataset.at(i),r[l][k],K));
			}	

			int bucket = cosine_g(gs,shuffle_idents[l],K);
			bucket = bucket % TableSize;

			assert(bucket >= 0);
			vector<Point>* myv = &(Hash_Tables[l][bucket]);
			dataset.at(i).g = shuffle_idents[l];
			myv->push_back(dataset.at(i));
		
		}
	}
	


	std::cout<<"----- Printing Hash -------"<<std::endl<<std::endl;

	for(int h = 0 ; h < L ; h++){
		std::cout <<"Hash Table "<< h << " : "<< std::endl;

		for(int b = 0 ; b < TableSize; b++){
			int bucket_num = 0;

			std::cout <<"\tBucket: "<< b << " : ";
			for(vector<Point>::iterator it = Hash_Tables[h][b].begin() ; 
				it != Hash_Tables[h][b].end() ;it++ ){

					bucket_num++;
			}
			std::cout<<bucket_num<<std::endl;
		}

	}


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
	output.open(output_path.c_str(),std::ofstream::out | std::ofstream::trunc);
	if (!output.is_open()){
		std::cout<<"Output path not correct"<<endl;
		exit(-2);
	}


	double  average_time = 0;
	double approx = 0;
	int found_items = 0;

	data_type max_ratio = -1;
	

	for(int query = 0 ; query < query_points; query++){

		output<<std::endl<<"--------------------"<<std::endl;
		output<<std::endl<<"--------------------"<<std::endl;
		output<<"(cosine)Query:"<< query<<std::endl;
		std::cout<<"(cosine)Query:"<< query<<std::endl;
		data_type min_dist = 10*R;
		const clock_t begin_time = clock();
		Point p = queries.at(query); 
		int items_found = 0;
		
		int nearest_id = -1;
		
		output<<"R-near neighbours:"<<std::endl;
		for(int l = 0 ; l < L ; l++){

			string gs;	
			for(int k = 0 ; k < K  ; k++){
				gs.push_back(cosine_h(p,r[l][k],K));
			}

			int bucket = cosine_g(gs,shuffle_idents[l],K);	
			bucket = bucket % TableSize;
			assert(bucket >= 0);
			for(vector<Point>::iterator it = (Hash_Tables[l][bucket]).begin() ; 
				it != (Hash_Tables[l][bucket]).end() ;it++ ){

				if(shuffle_idents[l] != it->g){
					continue;
				}

				data_type dist = CosineDistance(*it,p);				
				assert(dist >= 0);
				if(dist <= (1+e)*R){			

					items_found++;
					output<<"Item"<<it->get_id()<<std::endl;

					if(dist < min_dist){
						min_dist = dist;
						nearest_id = it->get_id();
					}

					if(items_found > 3*L)break;

				}
				
			}		
		}
		
		


		float tLSH =  (float( clock () - begin_time ) /  CLOCKS_PER_SEC);
		
		
		average_time += tLSH;	
		data_type true_min_dist=10*R;
		const clock_t begin_time_2 = clock();

		for(int i = 0 ; i < n; i++){
			data_type dist = CosineDistance(p,dataset.at(i));
			if(dist < true_min_dist)true_min_dist = dist;
		}
		
		output<<"Nearest neighbour: Item  ";
		if(nearest_id != -1){
			output<<nearest_id<<std::endl;
			output<<"distanceLSH :" << min_dist << std::endl;
		}
		else{
			output<<"not found"<<std::endl;
		}

		output<<"distanceTRUE :" << true_min_dist << std::endl;
		output<<"tLSH: "<<tLSH<<" sec"<<endl;
		output<<"tTrue:"<<(float( clock () - begin_time_2 ) /  CLOCKS_PER_SEC)<<" sec"<<endl;
		



		assert(min_dist >= true_min_dist);
		if(nearest_id != -1){
			data_type p_error = min_dist/true_min_dist;
			approx += p_error;
			if(p_error >= max_ratio)max_ratio = p_error;
			found_items++;
		}


		
	}

	output<<"*********************"<<std::endl;
	output<<"Average Query Time : "<<average_time/((float)query_points) << " seconds"<<std::endl;
	output<<"Mean LSH ratio: "<<approx/(float)found_items<<std::endl;
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
	int K = 4;
	int d = 0;
	int L = 5;
	int W = 350;
	int n;
	long long M;

	//Initialization
	if(argc == 1){
		input_path.assign("./input._mall");
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
		else if((s == "-L")and(argc > i+1)){
			L = atoi(argv[i+1]);
		}
		else if((s == "-o")and(argc > i+1)){
			string p = argv[i+1];
			output_path.assign(p);
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
			if( n == 0){//count d
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
		EuclideanLSH(n,L,K,d,W,query_path,dataset,output_path);
	else if(metric == "cosine")
		CosineLSH(n,dataset,K,L,d,query_path,output_path);
	else
		EuclideanLSH(n,L,K,d,W,query_path,dataset,output_path);

	return 0;



}