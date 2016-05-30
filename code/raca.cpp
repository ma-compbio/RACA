#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>
#include <list>
#include <iomanip>
#include <cmath>

using namespace std;

enum {SUCCESS, FAIL, CYCLE}; 
const double gMaxScore = 1.0;
const double gMinScore = 0.0000001;
const double gDiffScore = gMaxScore - gMinScore;
double gAlpha = 0.5;
int gIgnoreAdjs_woreads = 0;
int gMinPEReads = 100;	
double gMinRelCov = 0.5;

class Edge {
public:
	int bid1, bid2;
	int dir1, dir2;
	double score1, score2, weight;

	Edge(int _bid1, int _dir1, int _bid2, int _dir2)
			:bid1(_bid1), dir1(_dir1), bid2(_bid2), dir2(_dir2), 
			score1(0.0), score2(0.0), weight(0.0) {}
	
	bool operator<(const Edge& p) const {
		if (bid1 == p.bid1 && dir1 == p.dir1 && bid2 == p.bid2) return dir2 < p.dir2;
		if (bid1 == p.bid1 && dir1 == p.dir1) return bid2 < p.bid2;
		if (bid1 == p.bid1) return dir1 < p.dir1;
		return bid1 < p.bid1;
	}

	void reverse() {
		int tmp = bid1;
		bid1 = bid2;
		bid2 = tmp;

		tmp = dir1;
		dir1 = dir2;
		dir2 = tmp;
		if (dir1 == 1) dir1 = -1;
		else dir1 = 1;
		if (dir2 == 1) dir2 = -1;
		else dir2 = 1;
	}

	string toString(map<int, string>& map2name) {
		string bname1 = map2name[bid1];
		string bname2 = map2name[bid2];
		string strdir1 = "+";
		if (dir1 == -1) { strdir1 = "-"; }
		string strdir2 = "+";
		if (dir2 == -1) { strdir2 = "-"; }
		stringstream ss;
		ss << bname1 << " " << strdir1 << "\t" << bname2 << " " << strdir2 << "\t";
		ss << fixed << setprecision(10) << weight << "\t" << score1 << "\t" << score2;
		return ss.str(); 
	}
};

double compute_weight(int bid1, int dir1, int bid2, int dir2, map<Edge,double>& mapscore1, map<Edge,double>& mapscore2) 
{
	double adjscore = 0.0;
	double covscore = 0.0;
	map<Edge, double>::iterator bpiter1, bpiter2;
	map<Edge, int>::iterator iter;

	bpiter2 = mapscore2.find(Edge(bid1, dir1, bid2, dir2));
    if (bpiter2 != mapscore2.end()) covscore = bpiter2->second;

	if (gIgnoreAdjs_woreads == 1 && bpiter2 == mapscore2.end()) return 0.0;

	bpiter1 = mapscore1.find(Edge(bid1, dir1, bid2, dir2));
    if (bpiter1 != mapscore1.end()) adjscore = bpiter1->second;

	double weight = gAlpha*adjscore + (1-gAlpha)*covscore;	
	return weight;
}

void error (string msg, string file="") 
{
	cerr << msg << file << endl;
	exit(1);
}

bool cmp(const pair<Edge,double> &p1, const pair<Edge,double> &p2) 
{
	return p1.second > p2.second;
}

int insertEdge(list<Edge>& le, Edge& e, map<int,int>& mapUsed)
{
	Edge& fe = le.front();
	Edge& be = le.back();
	
	if (fe.bid1 == e.bid1 && fe.dir1 != e.dir1) {
		// check for a cycle
		if (be.bid2 == e.bid2 && be.dir2 != e.dir2) {
			return CYCLE; 
		}
		
		// e precedes fe with an opposite direction
		e.reverse();
		le.push_front(e);
		mapUsed[fe.bid1] = 1;
		mapUsed[-fe.bid1] = 1;
		if (e.dir1 == 1) mapUsed[-e.bid1] = 1;
		else mapUsed[e.bid1] = 1;
		return SUCCESS;	
	}
	if (fe.bid1 == e.bid2 && fe.dir1 == e.dir2) {
		
		// check for a cycle
		if (be.bid2 == e.bid1 && be.dir2 == e.dir1) return CYCLE; 

		// e precedes fe with the same direction
		le.push_front(e);
		mapUsed[fe.bid1] = 1;
		mapUsed[-fe.bid1] = 1;
		if (e.dir1 == 1) mapUsed[-e.bid1] = 1;
		else mapUsed[e.bid1] = 1;
		return SUCCESS;
	}

	if (be.bid2 == e.bid1 && be.dir2 == e.dir1) {
		// check for a cycle
		if (fe.bid1 == e.bid2 && fe.dir1 == e.dir2) return CYCLE; 
		
		// be precedes e with the same direction
		le.push_back(e);
		mapUsed[be.bid2] = 1;
		mapUsed[-be.bid2] = 1;
		if (e.dir2 == 1) mapUsed[e.bid2] = 1;
		else mapUsed[-e.bid2] = 1;
		return SUCCESS;
	}
	if (be.bid2 == e.bid2 && be.dir2 != e.dir2) {
		// check for a cycle
		if (fe.bid1 == e.bid1 && fe.dir1 != e.dir1) return CYCLE; 
		
		// be precedes e with an opposite direction
		e.reverse();
		le.push_back(e);
		mapUsed[be.bid2] = 1;
		mapUsed[-be.bid2] = 1;
		if (e.dir2 == 1) mapUsed[e.bid2] = 1;
		else mapUsed[-e.bid2] = 1;
		return SUCCESS;	
	}
	return FAIL;
}
	
void printLists(map<int,string>& mapBlockName, map<int, list<Edge> >& mapClasses)
{
	map<int, list<Edge> >::iterator citer;

	int clsnum = 1;	
	for(citer = mapClasses.begin(); citer != mapClasses.end(); citer++) {
		list<Edge>& le = citer->second;
		cout << ">" << clsnum << endl;
		clsnum++;
	
		list<Edge>::iterator liter;
		for (liter = le.begin(); liter != le.end(); liter++) {
			Edge& e = *liter;
			cout << e.toString(mapBlockName) << endl; 
		}	
		cout << endl;
	} // end of for
}

double computeAvg(map<Edge,double>& mapScores)
{
	map<Edge,double>::iterator iter;
	double sum = 0.0;
	int cnt = mapScores.size();
	for (iter = mapScores.begin(); iter != mapScores.end(); iter++) {
		sum += iter->second;
	}
	return (sum/cnt);
} 

double computeMedian(map<Edge,double>& mapScores, bool bIgnoreMin=false) 
{
    double midScore = 0.0;
    vector<pair<Edge,double> > vecScores(mapScores.begin(), mapScores.end());
    sort (vecScores.begin(), vecScores.end(), cmp);

	if (bIgnoreMin) {
		vector<pair<Edge,double> > vecScoresNew;
		for (int i = 0; i < vecScores.size(); i++) {
			pair<Edge,double> p = vecScores.at(i);
			Edge e = p.first;
			double score = p.second;
			if (score <= gMinScore) break;
			vecScoresNew.push_back(make_pair(e,score));
		}

		vecScores = vecScoresNew;	
	}

    int numEdges = vecScores.size();
    if (numEdges % 2 == 0) {
		int midindex2 = numEdges/2;
		int midindex1 = midindex2-1;
        pair<Edge,double> p = vecScores.at(midindex1);
        midScore = p.second;
        p = vecScores.at(midindex2);
        midScore += p.second;
		midScore /= 2;	
    } else {
        int midindex = (int)floor(numEdges/2);
        pair<Edge,double> p = vecScores.at(midindex);
        midScore = p.second;
    }

	return midScore;
}

void mergeLists(int clscnt, map<int, list<Edge> > &mapClasses) 
{
	map<int, list<Edge> >::iterator citer;
	list<Edge>::iterator liter;
	list<Edge>::reverse_iterator rliter;
	bool changed = true;
	while(changed) {	// repeat until no more changes
		changed = false;
		
		for (int i = 1; i <= clscnt; i++) {
			citer = mapClasses.find(i);
			if (citer == mapClasses.end()) continue;
			list<Edge>& le1 = citer->second;
		
			for (int j = 1; j <= clscnt; j++) {
				if (i == j) continue;

				citer = mapClasses.find(j);
				if (citer == mapClasses.end()) continue;
				list<Edge>& le2 = citer->second;
			
				Edge& e1front = le1.front();
				Edge& e1back = le1.back();
		
				Edge& e2front = le2.front();
				Edge& e2back = le2.back();

				if(e1front.bid1 == e2front.bid1 && e1front.dir1 != e2front.dir1) {
					// check cycle
					if (e1back.bid2 != e2back.bid2) {	
						for (liter = le2.begin(); liter != le2.end(); liter++) {
							Edge e2 = *liter;
							e2.reverse();
							le1.push_front(e2);
						} // end of for
						mapClasses.erase(j);
						changed = true;
					}
				} else if(e1front.bid1 == e2back.bid2 && e1front.dir1 == e2back.dir2) {
					// check cycle
					if (e1back.bid2 != e2front.bid1) {
						for (rliter = le2.rbegin(); rliter != le2.rend(); rliter++) {
							Edge e2 = *rliter;
							le1.push_front(e2);
						} // end of for
						mapClasses.erase(j);	
						changed = true;
					}	
				} else if(e1back.bid2 == e2front.bid1 && e1back.dir2 == e2front.dir1) {
					// check cycle
					if (e1front.bid1 != e2back.bid2) {
						for (liter = le2.begin(); liter != le2.end(); liter++) {
							Edge e2 = *liter;
							le1.push_back(e2);
						} // end of for
						mapClasses.erase(j);	
						changed = true;
					}	
				} else if(e1back.bid2 == e2back.bid2 && e1back.dir2 != e2back.dir2) {
					// check cycle
					if (e1front.bid1 != e2front.bid1) {
						for (rliter = le2.rbegin(); rliter != le2.rend(); rliter++) {
							Edge e2 = *rliter;
							e2.reverse();
							le1.push_back(e2);
						} // end of for
						mapClasses.erase(j);	
						changed = true;	
					}
				}	
			} // end of for j
		} // end of for i
	}
}

int main(int argc, char* argv[]) 
{
	gIgnoreAdjs_woreads = atoi(argv[1]);
	gAlpha = atof(argv[2]);
	char* fblist = argv[3];
	char* fcons = argv[4];
	char* flink = argv[5];	
	map<Edge, double>::iterator bpiter;

	if (gIgnoreAdjs_woreads == 1) cerr << "Ignore adjs. without read support? Yes" << endl;
	cerr << "Parameter alpha = " << gAlpha << endl;	
	cerr << "Block list = " << fblist << endl;
	cerr << "Conservation score file = " << fcons << endl;
	cerr << "Link score file = " << flink << endl;

	// read a list of blocks
	map<int, string> mapBlockName;
	map<string, int> mapBlockIndex;
	vector<string> vecBlocks;
	ifstream infile;
	infile.open(fblist);
	if (!infile) error ("\n[ERROR] Unable to open file: ", fblist); 
	int scfnum, start, end, blockid;
	string dir, rchr;
	int bindex = 1;
	while (infile >> scfnum >> start >> end >> dir >> blockid >> rchr) { 
		stringstream ss;
		ss << scfnum << ":" << start << "-" << end;
		string bname = ss.str();
		vecBlocks.push_back(bname);	
		mapBlockName[bindex] = bname;
		mapBlockIndex[bname] = bindex;
		bindex++;
	} // end of while	
	infile.close();

	// read adjacency scores
	map<Edge, double> mapAdjScores;
	infile.open(fcons);
	if (!infile) error ("\n[ERROR] Unable to open file: ", fcons); 
	string block1, dir1, block2, dir2;
	double adjscore;
	while (infile >> block1 >> dir1 >> block2 >> dir2 >> adjscore) {
		int bindex1 = mapBlockIndex[block1]; 
		int bindex2 = mapBlockIndex[block2];
		int intdir1 = 1;
		if (dir1.compare("-") == 0) intdir1 = -1; 
		int intdir2 = 1;
		if (dir2.compare("-") == 0) intdir2 = -1; 
		Edge bpf(bindex1, intdir1, bindex2, intdir2);
		Edge bpr(bindex2, -intdir2, bindex1, -intdir1);
		mapAdjScores[bpf] = adjscore;
		mapAdjScores[bpr] = adjscore;
	}
	infile.close();	

	// read coverage 
	map<Edge, double> mapCovScores;
	infile.open(flink);
	if (!infile) error ("\n[ERROR] Unable to open file: ", flink); 
	double covscore;
	while (infile >> block1 >> dir1 >> block2 >> dir2 >> covscore) {
		int bindex1 = mapBlockIndex[block1]; 
		int bindex2 = mapBlockIndex[block2];
		int intdir1 = 1;
		if (dir1.compare("-") == 0) intdir1 = -1; 
		int intdir2 = 1;
		if (dir2.compare("-") == 0) intdir2 = -1; 
		Edge bpf(bindex1, intdir1, bindex2, intdir2);
		Edge bpr(bindex2, -intdir2, bindex1, -intdir1);
		mapCovScores[bpf] = covscore;
		mapCovScores[bpr] = covscore;
	}
	infile.close();	
	
	// compute edge weights
	map<Edge, double> mapWeights;
	double weight = 0.0;
	int numblocks = vecBlocks.size();
	for (int i = 1; i <= numblocks; i++) {	// start from 1
		string bname1 = mapBlockName[i];
		for (int j = 1; j <= numblocks; j++) { // start from 1
			if (i == j) continue;
			string bname2 = mapBlockName[j];

			weight = compute_weight(i,1,j,1,mapAdjScores,mapCovScores);
			if (weight > 0.0) mapWeights[Edge(i,1,j,1)] = weight; 
			weight = compute_weight(i,1,j,-1,mapAdjScores,mapCovScores);
			if (weight > 0.0) mapWeights[Edge(i,1,j,-1)] = weight; 
			weight = compute_weight(i,-1,j,1,mapAdjScores,mapCovScores);
			if (weight > 0.0) mapWeights[Edge(i,-1,j,1)] = weight; 
			weight = compute_weight(i,-1,j,-1,mapAdjScores,mapCovScores);
			if (weight > 0.0) mapWeights[Edge(i,-1,j,-1)] = weight; 
		} // end of j
	} // end of for

	// greedy search based on edge weights
	vector<pair<Edge,double> > vecEdges(mapWeights.begin(), mapWeights.end()); 
	sort (vecEdges.begin(), vecEdges.end(), cmp);

	map<int, int> mapUsed;	// pos:head is connected, neg:tail is connected
	map<int, int>::iterator uiter, uiterex;
	map<int, list<Edge> > mapClasses;
	map<int, list<Edge> >::iterator citer;
	map<Edge, int>::iterator niter;
	map<Edge, double>::iterator biter;
	int clscnt = 0;
	for (int i = 0; i < vecEdges.size(); i++) {
		pair<Edge,double> p = vecEdges.at(i);
		Edge& e = p.first;
		e.weight = p.second;

		bpiter = mapAdjScores.find(e);
		if (bpiter != mapAdjScores.end()) e.score1 = mapAdjScores[e];
		bpiter = mapCovScores.find(e);
		if (bpiter != mapCovScores.end()) e.score2 = mapCovScores[e];
	
		if (e.weight < 0.1) continue;

		if (e.dir1 == 1) uiter = mapUsed.find(-e.bid1); 
		else uiter = mapUsed.find(e.bid1);  
		if (uiter != mapUsed.end()) {
			continue;				
		}
	
		if (e.dir2 == 1) uiter = mapUsed.find(e.bid2);
		else uiter = mapUsed.find(-e.bid2); 
		if (uiter != mapUsed.end()) {
			continue;		
		}		

		bool found = false;
		for (citer = mapClasses.begin(); citer != mapClasses.end(); citer++) {
			list<Edge>& le = citer->second;
			int res = insertEdge(le, e, mapUsed);
			if (res == SUCCESS || res == CYCLE) { 
				if (res == SUCCESS) mergeLists(clscnt, mapClasses);
				found = true;
				break;
			} // end of if	
		} // end of for 
	
		if (found == false) {
			list<Edge> le;
			le.push_back(e);
			mapClasses[++clscnt] = le;
		
			if (e.dir1 == 1) mapUsed[-e.bid1] = 1; 
			else mapUsed[e.bid1] = 1; 
			if (e.dir2 == 1) mapUsed[e.bid2] = 1; 
			else mapUsed[-e.bid2] = 1; 
		}	
	}

	mergeLists(clscnt, mapClasses);
	printLists(mapBlockName, mapClasses);

	return 0;
}
