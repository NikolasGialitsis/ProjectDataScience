#ifndef LSH_CLASSES_H
#define LSH_CLASSES_H


typedef double data_type ;
typedef std::string bitstring;

using namespace std;
int id = 0;
int space = 0;

class Point{
public:

	std::vector<int> g;
	Point(){
		cout << "Construct a default point"<<endl;
		this->item_id = ++id;
		this->fields.clear();		
		(this->g).clear();

	}
	Point(string s){

		/*int loc = s.find("\t");
		string temp = s.substr(0,loc);
		this->item_id = atoi(temp.c_str());
		this->fields.assign(s);
		*/
		this->item_id = id++;
		space += sizeof(int);
		space += s.size();
	}

		int get_id(){
			return this->item_id;
		};
		
		string get_all_fields(){
			return this->fields;
		};

private:
	int item_id;
	string fields;
	
};



void EuclideanLSH(int n,int L,int K,int d,int W,string query_path,std::vector<Point>& dataset,string& output_path);
double EuclideanDistance(Point& p1,Point& p2);
int euclidean_h(Point& point,int w,int d,std::vector<data_type>& v, data_type t);
void euclidean_g(std::vector<int>& ,std::vector<int>& h,std::vector<int>& shuffle_idents,int K);
int euclidean_phi(std::vector<int>& r, std::vector<int>& g,int K,int TableSize);

void CosineLSH(int n,std::vector<Point>& dataset,int K,int L,int d,string& query_path,string& output_path);
char cosine_h(Point& point,vector<data_type>& r,int K);
void cosine_g(bitstring& gv,int h[],std::vector<int>& shuffle,int K);
double CosineDistance(Point& p1,Point& p2);

void CosineCube(int n,std::vector<Point>& dataset,int K,int d,string& query_path,string& output_path);
void EuclideanCube(int n,int K,int d,int W,string query_path,std::vector<Point>& dataset,string& output_path);
int euclidean_h(Point& point,int w,int d,std::vector<data_type>& v,data_type t);
char cosine_h(Point& point,vector<data_type>& r,int K);
string euclidean_f(std::vector<int> h);



int hammingDistance(string& s1,string& s2);
int BinaryToDecimal(string& h);
int modulo(int x,int y);
int long_modulo(int x,long long y);
unsigned b2d(unsigned num);
int BinaryToDecimal(string& h);

#endif
