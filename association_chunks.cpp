#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>

#include <chrono>

using namespace std;


/*

One area common to EROS and MACHO. Associating stars there. (center "ori")

*/

class Point{
public:
	Point(double, double);
	Point(string, string);
	double angular_distance(Point&);
	void projection(double alpha0=1.4, double delta0=-1.2);
	//~Point(){cout << _identifier;}
// variables
	double _alpha_rad, _delta_rad;
	double _alpha_deg, _delta_deg;
	double _a, _d;
	string _ra, _dec;
	string _identifier;
	string _red_chunk;
	bool _mod = false;
	vector<Point*> _nearest;
};

typedef vector<Point*> pointSet;
typedef map<int, pointSet> pointSetMap;
typedef map<string, pointSet> chunkMap;

vector<string> split(string str, char delimiter = ' ');
double rad_to_deg(double);
double deg_to_rad(double);
double sign(double);
double right_ascension_to_rad(vector<string>);
double declination_to_rad(vector<string>);
vector<double> shortest_distances(vector<Point>&, Point&, int vector_size=3);
void closest_points(vector< pair<double, Point> >& holder, vector<Point>& pts_array, Point& origin, double vector_size=3);
bool pair_point_compare(pair<double, Point>& p1, pair<double, Point>& p2);
void load_MACHO_stars(vector<Point>& stars, int field_id);
void load_EROS_stars(vector<Point>& stars, string field);
void load_OGLE3_stars(vector<Point>& stars, int field);
void find_stars_in_area(vector<Point>& selected_ones, vector<Point>& stars, Point& origin, double alpha=0.0035, double delta=0.0013, double mininum_separation_distance=1.0e-5, int identifier_range=1);
void rec_place(Point& origin, Point& to_insert, int i=0);
void find_nearest_stars(vector<Point>& stars1, vector<Point>& stars2, int vector_size=3);
void find_nearest_stars(vector<Point*>& stars1, vector<Point*>& stars2, int vector_size=3);
double projection(Point& point, double alpha0=1.4);
int which_cell(Point& point, Point& origin, int nb_cols, int nb_rows, double width, double height);
void divide_in_cells(pointSetMap& cells, vector<Point>& stars, Point& origin, int nb_cols, int nb_rows, double width, double height);
void divide_in_cells(pointSetMap& cells, pointSet& stars, Point& origin, int nb_cols, int nb_rows, double width, double height);
void associate(pointSetMap& cells1, pointSetMap& cells2, int nb_cols, int nb_rows);

int correct_position(pointSet& stars);
void chunkify(vector<Point>& stars, chunkMap& macho_chunks);

double mean(vector<double>& vect);
double std_deviation(vector<double>& vect);

int main(){
	cout << "Starting" << endl;
	chrono::time_point<chrono::system_clock> start, loading, end;
	start = chrono::system_clock::now();

	// int nb_cols = 483, nb_rows=480;
	int nb_cols = 966, nb_rows=960;
	Point origin(0.0, 0.0);
	origin._a = 84.6990*M_PI/180.-1.4;
	origin._d = -64.9931*M_PI/180.+1.2;

	vector<Point> eros_stars, macho_stars;
	chunkMap macho_chunks;
	load_EROS_stars(eros_stars, "060");
	load_MACHO_stars(macho_stars, 15);
	chunkify(macho_stars, macho_chunks);
	pointSetMap eros_cells, macho_cells;

	loading = chrono::system_clock::now();

	divide_in_cells(eros_cells, eros_stars, origin, nb_cols, nb_rows, 0.140485, 0.139755);
	cout << "Looping through chunks." << endl;
	bool do_continue;
	int status, counter;
	for(chunkMap::iterator it=macho_chunks.begin(); it!=macho_chunks.end(); it++)
	{
		do_continue = true;
		status = 0;
		counter = 0;
		while(do_continue)
		{
			counter++;
			//divide chunks stars in cells
			macho_cells.clear();
			divide_in_cells(macho_cells, it->second, origin, nb_cols, nb_rows, 0.140485, 0.139755);
			//associate with eros stars
			associate(eros_cells, macho_cells, nb_cols, nb_rows);
			//correct position of macho stars
			status = correct_position(it->second);
			cout << status << endl;
			if(status or counter>20) do_continue=false;
		}
		cout << it->first << endl;
	}

	//generate closest stars from corrected macho_stars
	macho_cells.clear();
	divide_in_cells(macho_cells, macho_stars, origin, nb_cols, nb_rows, 0.140485, 0.139755);
	associate(eros_cells, macho_cells, nb_cols, nb_rows);

	cout << eros_stars.size() << endl;

	cout << "Selection closer-closer" << endl;

	ofstream out("./temp.txt");
	ofstream out1("./idlist.txt");
	for(vector<Point>::iterator it = eros_stars.begin(); it!=eros_stars.end(); it++)
	{
		if(not it->_nearest.empty() and not it->_nearest[0]->_nearest.empty())
		{	
			if(it->_identifier == it->_nearest[0]->_nearest[0]->_identifier)// and it->_nearest[1]->_nearest[0]->_identifier == it->_nearest[0]->_nearest[1]->_identifier)
			{
				out << setprecision(9) << it->_alpha_rad << " " << it->_delta_rad << " " << it->_nearest[0]->_alpha_rad << " "  << it->_nearest[0]->_delta_rad << " " << it->_nearest[1]->_alpha_rad << " " << it->_nearest[1]->_delta_rad << endl;
				out1 << it->_identifier << " " << it->_nearest[0]->_identifier << endl;
			}
		}	
	}
	out.close();
	out1.close();

	end = chrono::system_clock::now();

	cout << "Finished." << endl;
	cout << "Loading time : " << chrono::duration_cast<chrono::milliseconds>(loading-start).count() << " ms" << endl;
	cout << "Computing time : " << chrono::duration_cast<chrono::milliseconds>(end-loading).count() << " ms" << endl;
	return 0;
}

vector<string> split(string str, char delimiter)
	/**
	Split a string into words seperated by a delimiter

	@param str: 		The string to be split
	@param delimiter: 	The character where the string will be split

	@return: A vector containing all the resulting pieces
	*/
{
	stringstream strh(str);
	string temp;
	vector<string> v;
	while(getline(strh, temp, delimiter))
	{
		if(temp!="") v.push_back(temp);
	}
	return v;
}

double rad_to_deg(double x)
	/**
	Convert an angle from radians to degrees.
	*/
{
	return x*180./M_PI;
}

double deg_to_rad(double x)
	/**
	Convert an angle from degrees to radians.
	*/
{
	return x*M_PI/180.;
}

double sign(double x)
{
	if(x != 0)
	{
		return x/abs(x);
	}
	return 1;
}

double right_ascension_to_rad(string ra_str)
	/**
	Convert a right ascension in format h:min:sec into alpha in radians

	@param ra_str: right ascension "h:min:sec"
	@return: alpha in radians
	*/
{
	vector<string> ra_str_vec;
	ra_str_vec= split(ra_str, ':');
	try{
		if(ra_str_vec.size()!=3){
			throw string("Bad dimension for right ascension!");
		}
		else{
			vector<double> ra;
			for(int i=0;i<3;i++){
				ra.push_back(stod(ra_str_vec[i]));
			}
			return 2.*M_PI/24.*(ra[0] + (ra[1] + ra[2]/60.)/60.); 
		}
	}
	catch(string& err){
		cerr << err << endl;
	}
	return 0;
}

double declination_to_rad(string dec_str)
	/**
	Convert a declination in format deg:min:sec into delta in radians

	@param ra_str: right ascension "deg:min:sec"
	@return: delta in radians
	*/
{
	vector<string> dec_str_vec;
	dec_str_vec= split(dec_str, ':');
	try{
		if(dec_str_vec.size()!=3){
			throw string("Bad dimension for declination!");
		}
		else{
			vector<double> dec;
			for(int i=0;i<3;i++){
				dec.push_back(stod(dec_str_vec[i]));
			}
			return 2.*M_PI/360.*(abs(dec[0]) + (dec[1] + dec[2]/60.)/60.) * sign(dec[0]);
		}
	}
	catch(string& err){
		cerr << err << endl;
	}
	return 0;
}

void closest_points(vector< pair<double,Point> >& holder, vector<Point>& pts_array, Point& origin, double vector_size)
{
	/*
	Find the "vector_size" closest points to "origin" in list "pts_array", and put them into "holder", with the associated distance.
	*/
	holder.clear();
	double distance;
	for(int i=0;i<pts_array.size();i++)
	{
		distance = origin.angular_distance(pts_array[i]);
		if(holder.size()<vector_size or distance < holder[2].first)
		{
			holder.push_back(make_pair(distance, pts_array[i]));
			if(holder.size()==vector_size+1)
			{
				sort(holder.begin(), holder.end(), pair_point_compare);
				holder.pop_back();
			}
		}
	}
}

vector<double> shortest_distances(vector<Point>& stars_centers, Point& centre, int vector_size)
	/*
	Find the shortest distances between "centre" and all the points in "star_centers".
	*/
{
	vector<double> distances;
	for(int i=0; i<stars_centers.size(); i++)
	{
		double curr_dist = centre.angular_distance(stars_centers[i]);
		if(distances.size()<vector_size or curr_dist<distances[vector_size-1])
			distances.push_back(curr_dist);
			if(distances.size()<vector_size)
			{
				sort(distances.begin(), distances.end());
				distances.pop_back();
			}
	}
	return distances;
}

bool pair_point_compare(pair<double, Point>& p1, pair<double, Point>& p2)
{
	return p1.first < p2.first;
}

void load_MACHO_stars(vector<Point>& stars, int field_id)
{
	/*
	Load stars coordinates and identifier of a MACHO field F_"fieldid" in "stars".
	*/
	stringstream filename;
	string line;
	vector<string> splitline;
	filename << "/Users/tristanblaineau/Documents/MACHO/Positions/DumpStar_" << field_id << ".txt";
	cout << "Loading " << filename.str() << endl;
	ifstream f(filename.str());
	if(f.is_open())
	{
		while(getline(f, line))
		{
			splitline = split(line, ';');
			if(true)//splitline[6]+splitline[7]=="W18")
			{
				Point p(splitline[3], splitline[4]);
				p._identifier = splitline[0]+":"+splitline[1]+":"+splitline[2];
				p._red_chunk = splitline[6]+splitline[7];
				stars.push_back(p);
			}
		}
	}
	else cout << "Could not load file." << endl;
	f.close();
}

void chunkify(vector<Point>& stars, chunkMap& macho_chunks)
{
	for(vector<Point>::iterator it=stars.begin(); it!=stars.end(); it++)
	{
		macho_chunks.insert(pair<string, pointSet>(it->_red_chunk, pointSet()));
		macho_chunks[it->_red_chunk].push_back(&*it);
	}
}

void load_EROS_stars(vector<Point>& stars, string field)
{
	/*
	Load stars coordinates and identifier of a EROS subfield.
	*/

	stringstream filename;
	string path = "/Volumes/DisqueSauvegarde/EROS/lm/";
	string ccds[] = {"k", "l", "m", "n"};
	vector<string> splitline;
	string line;
	for(int i=0; i<=7; i++)
	{
		for(int j=0; j<=3; j++)
		{
			filename.str("");
			filename << path << "lm" << field << "/lm" << field << to_string(i) << "/lm" << field << to_string(i) << ccds[j] << "/lm" << field << to_string(i) << ccds[j] << ".cat";
			cout << "Loading " << filename.str() << endl;
			ifstream f(filename.str());
			if(f.is_open())
			{
				getline(f, line);
				while(getline(f, line))
				{
					splitline = split(line, ' ');
					Point p(stod(splitline[1]), stod(splitline[2]));
					p._identifier = splitline[0];
					stars.push_back(p);
				}
			}
			else cout << "Could not load file." << endl;
			f.close();
		}
	}
}

void load_OGLE3_stars(vector<Point>& stars, int field)
{
	stringstream filename, identifier;
	string path = "/Volumes/DisqueSauvegarde/OGLE-III/lmc/";
	vector<string> splitline;;
	string line;
	for(int i=1; i<=8; i++)
	{
		filename.str("");
		filename << path << "lmc" << to_string(field) << "." << to_string(i) << ".map";
		cout << "Loading " << filename.str() << endl;
		ifstream f(filename.str());
		if(f.is_open())
		{
			while(getline(f, line))
			{
				identifier.str("");
				splitline = split(line, ' ');
				Point p(splitline[1], splitline[2]);
				identifier << to_string(field) << ":" << to_string(i) << ":" << splitline[0];
				p._identifier = identifier.str();
				stars.push_back(p);
			}
			f.close();
		}
		else cout << "Could not load file !" << endl;
	}
}

void find_stars_in_area(vector<Point>& selected_ones, vector<Point>& stars, Point& origin, double alpha, double delta, double mininum_separation_distance, int identifier_range)
{
	/*
	Returns in "selected_ones" all points from "stars" that are inside the rectangle of center "origin" and angular width and length "alpha","delta".

	selected_points: 				stocks the stars found in the area
	stars: 							scanned stars
	origin: 						the center of the zone we want
	alpha:
	delta: 							angular size of the area in radians
	minimum_separation_distance: 	minimum distance under which two stars from two different fields are considered to be the same
	identifier range: 				length of the characters that are needed to identify the field in the identifier of the star (ex: for MACHO, 1 "6:2344:12", "82:9087:789" -> "6:", "82" )
	*/
	double alpha_min = origin._alpha_rad - alpha/2;
	double alpha_max = origin._alpha_rad + alpha/2;
	double delta_min = origin._delta_rad - delta/2;
	double delta_max = origin._delta_rad + delta/2;
	double distance;
	bool can_add;
	for(int i=0; i<stars.size(); i++)
	{
		if(stars[i]._alpha_rad < alpha_max and 
			stars[i]._alpha_rad > alpha_min and
			stars[i]._delta_rad < delta_max and
			stars[i]._delta_rad > delta_min)
		{
			can_add = true;
			for (int j=0; j<selected_ones.size();j++)
			{
				// Checking if there are other stars of an overlapping field that could be the same one recorded twice.
				distance = selected_ones[j].angular_distance(stars[i]);
				if(distance<mininum_separation_distance and selected_ones[j]._identifier.substr(0,identifier_range)!=stars[i]._identifier.substr(0,identifier_range))
				{
					can_add = false;
					break;
				}
			}
			if(can_add) selected_ones.push_back(stars[i]);
		}
	}
}

void find_nearest_stars(vector<Point>& stars1, vector<Point>& stars2, int vector_size)
{
	/*
	For all stars in "stars1", find the "vector_size" nearest stars in "stars2"
	
	stars1 is the small field, stars2 is the big field, only stars1 are modified.
	*/

	double distance;
	for(int i=0; i<stars1.size(); i++)
	{
		for(int j=0; j<stars2.size(); j++)
		{
			if(stars1[i]._nearest.size()==0) stars1[i]._nearest.push_back(&stars2[j]);
			else
			{
				rec_place(stars1[i], stars2[j]);
				if(stars1[i]._nearest.size()>vector_size) stars1[i]._nearest.pop_back();
			}
		}
	}
}

void find_nearest_stars(vector<Point*>& stars1, vector<Point*>& stars2, int vector_size)
{
	/*
	For all stars in "stars1", find the "vector_size" nearest stars in "stars2"
	
	stars1 is the small field, stars2 is the big field, only stars1 are modified.
	*/

	double distance;
	for(int i=0; i<stars1.size(); i++)
	{
		for(int j=0; j<stars2.size(); j++)
		{
			if(stars1[i]->_nearest.size()==0) stars1[i]->_nearest.push_back(stars2[j]);
			else
			{
				rec_place(*stars1[i], *stars2[j]);
				if(stars1[i]->_nearest.size()>vector_size) stars1[i]->_nearest.pop_back();
			}
		}
	}
}

void rec_place(Point& origin, Point& to_insert, int i)
{
	int t = origin._nearest.size();
	if(i!=t and origin._nearest[t-i-1]->angular_distance(origin) > origin.angular_distance(to_insert))
	{
		rec_place(origin, to_insert, i+1);
	}
	else origin._nearest.insert(origin._nearest.end()-i, &to_insert);
}

int which_cell(Point& point, Point& origin, int nb_cols, int nb_rows, double width, double height)
{
	/*
	Associate a cell number to "point". The cells are obtained by dividing a rectangle of upper left corner "origin",
	width and height "width" and "height" in projected coordinates in "nb_cols" columns and "nb_rows" rows. The cells
	are numbered from left to right, up to down, such as n = col + row * nb_cols with col the indice of the column, and
	row the indice of the row, both starting from 0.

	point: 		point to choose a cell
	origin:		upper left corner of the rectangle
	nb_cols: 	number of columns
	nb_row: 	number of rows
	width:		width of the rectangle, in radians, in the projected coordinates
	height:		height of the rectangle, in radians, in the projected coordinates

	return:		the number of the cell in which the point is situated.
	*/
	point.projection();
	double ntgprt;

	modf((origin._d - point._d)/height*nb_rows, &ntgprt);
	int row = ntgprt;
	if(row < 0 or row >= nb_rows)
	{
		// cout << "Invalid row : " << row << endl;
		return -1;
	}
	modf((origin._a - point._a)/width*nb_cols, &ntgprt);
	int col = ntgprt;
	if(col < 0 or col >= nb_cols)
	{
		// cout << "Invalid col : " << col << endl;
		return -1;
	}
	return col+row*nb_cols;
}

void divide_in_cells(pointSetMap& cells, vector<Point>& stars, Point& origin, int nb_cols, int nb_rows, double width, double height)
{
	/*
	Associate a cell number to each star in "stars" 
	*/
	int cell_nb;
	for(vector<Point>::iterator it = stars.begin(); it != stars.end(); it++)
	{
		cell_nb = which_cell(*it, origin, nb_cols, nb_rows, width, height);
		cells.insert(pair<int, pointSet>(cell_nb, pointSet()));
		cells[cell_nb].push_back(&*it);
	}
	cells.erase(-1); //Remove out of bounds stars (which_cells returned -1)
}

void divide_in_cells(pointSetMap& cells, vector<Point*>& stars, Point& origin, int nb_cols, int nb_rows, double width, double height)
{
	/*
	Associate a cell number to each star in "stars" 
	*/
	int cell_nb;
	for(vector<Point*>::iterator it = stars.begin(); it != stars.end(); it++)
	{
		cell_nb = which_cell(**it, origin, nb_cols, nb_rows, width, height);
		cells.insert(pair<int, pointSet>(cell_nb, pointSet()));
		cells[cell_nb].push_back(&**it);
	}
	cells.erase(-1);
}

int correct_position(pointSet& macho_stars)
	/*

	return: 0 if radius too wide
			1 if centered
			2 if too few stars
	*/
{
	vector<double> alpha_diff, delta_diff, alpha_diff2, delta_diff2;
	double mean_alpha, mean_delta, corr_alpha, corr_delta, radius;
	double mean_alpha2, mean_delta2, corr_alpha2, corr_delta2, radius2;
	int counter;
	if(macho_stars.size()>100)
	{
		mean_alpha = 0.;
		mean_delta = 0.;
		mean_alpha2 = 0.;
		mean_delta2 = 0.;
		counter = 0;
		for(int i=0; i<macho_stars.size(); i++)
		{
			if(macho_stars[i]->_nearest.size()>1 and macho_stars[i]->_identifier == macho_stars[i]->_nearest[0]->_nearest[0]->_identifier)
			{
				alpha_diff.push_back(macho_stars[i]->_alpha_rad - macho_stars[i]->_nearest[0]->_alpha_rad);
				delta_diff.push_back(macho_stars[i]->_delta_rad - macho_stars[i]->_nearest[0]->_delta_rad);

				alpha_diff2.push_back(macho_stars[i]->_alpha_rad - macho_stars[i]->_nearest[1]->_alpha_rad);
				delta_diff2.push_back(macho_stars[i]->_delta_rad - macho_stars[i]->_nearest[1]->_delta_rad);
			}
		}
		corr_alpha = mean(alpha_diff);
		corr_delta = mean(delta_diff);
		corr_alpha2 = mean(alpha_diff2);
		corr_delta2 = mean(delta_diff2);

		radius = sqrt(corr_alpha*corr_alpha+corr_delta*corr_delta);
		radius2 = sqrt(corr_alpha2*corr_alpha2+corr_delta2*corr_delta2);
		if(std_deviation(alpha_diff) > 1e-7 and std_deviation(delta_diff) > 1e-7)
		{
			corr_alpha = corr_alpha2;
			corr_delta = corr_delta2;
		}

		for(int i=0; i<macho_stars.size(); i++)
		{
			macho_stars[i]->_alpha_rad -= corr_alpha;
			macho_stars[i]->_delta_rad -= corr_delta;
			macho_stars[i]->_mod = true;
		}
		
		cout << radius*180*3600/M_PI << " " << radius2*180*3600/M_PI << endl;
		cout << std_deviation(alpha_diff)*180*3600/M_PI << " " << std_deviation(delta_diff)*180*3600/M_PI << " : " << std_deviation(alpha_diff2)*180*3600/M_PI << " " << std_deviation(delta_diff2)*180*3600/M_PI << endl;
		if(radius < 1e-7) return 1;
		return 0;
	}
	else 
		cout << "Too few stars in chunk to correct" << endl;
		return 2;
}

void associate(pointSetMap& cells1, pointSetMap& cells2, int nb_cols, int nb_rows)
{
	for(pointSetMap::iterator it=cells2.begin(); it!=cells2.end(); it++){
		for(int i=0; i<it->second.size();i++){
			it->second[i]->_nearest.clear();
		}
	}
	for(pointSetMap::iterator it=cells1.begin(); it!=cells1.end(); it++){
		for(int i=0; i<it->second.size();i++){
			it->second[i]->_nearest.clear();
		}
	}

	int index, bar_width=30;
	pointSet holder;
	double i=0., max=cells1.size(), prec=100, now=100;
	for(pointSetMap::iterator it = cells1.begin(); it != cells1.end(); it++)
	{
		// Loading bar
		modf(100*(max-i)/max, &now);
		cout << "[";
		for(int j=0; j<bar_width; j++)
		{
			if(j*100./bar_width < 100-prec) cout << "=";
			else cout << " ";
		}
		if(now != prec)
		{
			// cout << now << " " << i << endl;
			prec = now;	
		}
		cout << "] " << (100-now) <<" %\r";
		cout.flush();
		i++;

		index = it->first;
		// cout << index << endl;
		pointSetMap::iterator mcc = cells2.find(index);
		if(mcc != cells2.end()) holder = mcc->second;

		mcc = cells2.find(index-1);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells2.find(index+1);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells2.find(index+nb_cols-1);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells2.find(index+nb_cols);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells2.find(index+nb_cols+1);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells2.find(index-nb_cols-1);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells2.find(index-nb_cols);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells2.find(index-nb_cols+1);
		if(mcc != cells2.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		if(not holder.empty()){
			find_nearest_stars(it->second, holder);
			holder.clear();
		}
	}
	cout << endl;


	i=0.;
	max=cells2.size();
	prec=100;
	now=100;
	for(pointSetMap::iterator it = cells2.begin(); it != cells2.end(); it++)
	{

		// Loading bar
		modf(100*(max-i)/max, &now);
		cout << "[";
		for(int j=0; j<bar_width; j++)
		{
			if(j*100./bar_width < 100-prec) cout << "=";
			else cout << " ";
		}
		if(now != prec)
		{
			// cout << now << " " << i << endl;
			prec = now;	
		}
		cout << "] " << (100-now) <<" %\r";
		cout.flush();
		i++;


		index = it->first;
		// cout << index << endl;
		pointSetMap::iterator mcc = cells1.find(index);
		if(mcc != cells1.end()) holder = mcc->second;

		mcc = cells1.find(index-1);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells1.find(index+1);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells1.find(index+nb_cols-1);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells1.find(index+nb_cols);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells1.find(index+nb_cols+1);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells1.find(index-nb_cols-1);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells1.find(index-nb_cols);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		mcc = cells1.find(index-nb_cols+1);
		if(mcc != cells1.end()) holder.insert(holder.end(), mcc->second.begin(), mcc->second.end());

		if(not holder.empty()){
			find_nearest_stars(it->second, holder);
			holder.clear();
		}
	}
	cout << endl;
}

double mean(vector<double>& vect)
{
	double mean=0.;
	for(int i=0; i<vect.size(); i++)
	{
		mean+=vect[i];
	}
	return mean/vect.size();
}

double std_deviation(vector<double>& vect)
{
	double sd, mean1;
	sd=0;
	mean1 = mean(vect);
	for(int i=0; i<vect.size(); i++)
	{
		sd+=pow(vect[i]-mean1, 2);
	}
	return sd;
}

Point::Point(double adeg, double ddeg)
{
	_alpha_deg = adeg;
	_delta_deg = ddeg;
	_alpha_rad = deg_to_rad(adeg);
	_delta_rad = deg_to_rad(ddeg);
}

Point::Point(string raraw, string decraw)
{
	_ra = raraw;
	_dec = decraw;
	_alpha_rad = right_ascension_to_rad(raraw);
	_delta_rad = declination_to_rad(decraw);
}

double Point::angular_distance(Point& p)
	/*
	Compute angular distance in radians between current point and point p.
	*/
{
	double alpha_1, alpha_2, delta_1, delta_2;
	alpha_1 = _alpha_rad;
	delta_1 = _delta_rad;
	alpha_2 = p._alpha_rad;
	delta_2 = p._delta_rad;
	if(alpha_1==alpha_2 and delta_1==delta_2) return 0.;
	return acos(sin(delta_1)*sin(delta_2) + cos(delta_1)*cos(delta_2)*cos(alpha_1 - alpha_2));
}

void Point::projection(double alpha0, double delta0)
	/*
	Compute the projected coordinates (sphere to flat) with origin coordinates alpha0, delta0	
	*/
{
	double ndelta, nalpha, a;
	double deltag = M_PI/2.+delta0, alphag = alpha0;
	ndelta = asin(sin(deltag)*sin(_delta_rad) + cos(deltag)*cos(_delta_rad)*cos(_alpha_rad - alphag));
	nalpha = (M_PI - acos((sin(_delta_rad)-sin(deltag)*sin(ndelta))/(cos(deltag)*cos(ndelta))))*sign(_alpha_rad - alphag);
	a = acos(cos(nalpha)*cos(ndelta)*cos(ndelta)+sin(ndelta)*sin(ndelta))*sign(nalpha);
	_a = a;
	_d = ndelta;
}