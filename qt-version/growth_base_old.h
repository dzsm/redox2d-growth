// ///////////////////////////// MIT License //////////////////////////////////// //
//                                                                                //
// Copyright (c) 2010 David Zsolt Manrique                                        //
//                    david.zsolt.manrique@gmail.com                              //
//                                                                                //
// Permission is hereby granted, free of charge, to any person obtaining a copy   //
// of this software and associated documentation files (the "Software"), to deal  //
// in the Software without restriction, including without limitation the rights   //
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      //
// copies of the Software, and to permit persons to whom the Software is          //
// furnished to do so, subject to the following conditions:                       //
//                                                                                //
// The above copyright notice and this permission notice shall be included in     //
// all copies or substantial portions of the Software.                            //
//                                                                                //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     //
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         //
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  //
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN      //
// THE SOFTWARE.                                                                  //
//                                                                                //
// ////////////////////////////////////////////////////////////////////////////// //

#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <list>
#include <set>
#include <cstdarg>
#include <tuple>
#include <stdexcept>
#include <random>

#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <Eigen/Sparse>

//#include <cmdl_withorder.h>
//#include <ioutil.h>

//std::default_random_engine generator;
//std::uniform_real_distribution<double> distribution(0.0,1.0);

inline std::vector <std::string> &split(const std::string &s, char delim, std::vector <std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

inline std::vector <std::string> split(const std::string &s, char delim) {
    std::vector <std::string> elems;
    split(s, delim, elems);
    return elems;
}

//
//     2--1
//    / \/ \
//   3--12--0
//    \ /\ /
//     4--5
//

inline int index(int i, int j, int nw, int nh) {
    i = (i % nw + nw) % nw;
    j = (j % nh + nh) % nh;
    return j * nw + i;
}

inline void generateGrid(Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                         Eigen::Matrix<double, Eigen::Dynamic, 2> &co,
                         Eigen::Matrix<double, 6, 2> &an, int nw, int nh) {

    nb.resize(nw * nh, Eigen::NoChange);
    co.resize(nw * nh, Eigen::NoChange);

    const double sqrt3p2 = std::sqrt(3.0) / 2.0;

    an(0, 0) = 1.0;
    an(0, 1) = 0.0;
    an(1, 0) = 0.5;
    an(1, 1) = sqrt3p2;
    an(2, 0) = -0.5;
    an(2, 1) = sqrt3p2;
    an(3, 0) = -1.0;
    an(3, 1) = 0.0;
    an(4, 0) = -0.5;
    an(4, 1) = -sqrt3p2;
    an(5, 0) = 0.5;
    an(5, 1) = -sqrt3p2;

    for (int j = 0; j < nh; j++)
        for (int i = 0; i < nw; i++) {
            int k = j * nw + i;

            nb(k, 12) = k;

            if (j % 2 == 0) {
                nb(k, 0) = index(i + 1, j, nw, nh);
                nb(k, 1) = index(i, j + 1, nw, nh);
                nb(k, 2) = index(i - 1, j + 1, nw, nh);
                nb(k, 3) = index(i - 1, j, nw, nh);
                nb(k, 4) = index(i - 1, j - 1, nw, nh);
                nb(k, 5) = index(i, j - 1, nw, nh);

                co(k, 0) = i * 1.0;
                co(k, 1) = j * sqrt3p2;

            } else {
                nb(k, 0) = index(i + 1, j, nw, nh);
                nb(k, 1) = index(i + 1, j + 1, nw, nh);
                nb(k, 2) = index(i, j + 1, nw, nh);
                nb(k, 3) = index(i - 1, j, nw, nh);
                nb(k, 4) = index(i, j - 1, nw, nh);
                nb(k, 5) = index(i + 1, j - 1, nw, nh);

                co(k, 0) = i * 1.0 + 0.5;
                co(k, 1) = j * sqrt3p2;
            }
        }
    for (int k = 0; k < nb.rows(); k++) {
        nb(k, 6) = nb(nb(k, 0), 1);
        nb(k, 7) = nb(nb(k, 1), 2);
        nb(k, 8) = nb(nb(k, 2), 3);
        nb(k, 9) = nb(nb(k, 3), 4);
        nb(k, 10) = nb(nb(k, 4), 5);
        nb(k, 11) = nb(nb(k, 5), 0);
    }

}

inline void loadLabels(int &nw, int &nh, std::vector <std::string> &labels, std::istream &is, const std::string &key) {
    std::string findkey;
    while (findkey != key) is >> findkey;

    is >> nw >> nh;

    labels.resize(nw * nh);

    for (int j = nh - 1; j >= 0; j--)
        for (int i = 0; i < nw; i++) {
            int k = j * nw + i;
            is >> labels[k];
        }

}

template<typename T>
void sparseMatrixSolver(Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
                        const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                        const Eigen::Matrix<T, Eigen::Dynamic, 6> &e,
                        const Eigen::Matrix<T, Eigen::Dynamic, 1> &b) {

    Eigen::SparseMatrix <T> poisson(nb.rows(), nb.rows());         // default is column major
    poisson.reserve(Eigen::VectorXi::Constant(nb.rows(), 7));


    for (int k = 0; k < nb.rows(); k++) {
        for (int i = 0; i < 6; i++)
            poisson.insert(k, nb(k, i)) = e(k, i);
        poisson.insert(k, k) = -1.0;
    }
    poisson.makeCompressed();

    Eigen::BiCGSTAB <Eigen::SparseMatrix<T>> solver(poisson);

    u = solver.solveWithGuess(-b, u);

}


struct Label {
    int nw, nh;
    std::vector <std::string> label;

    inline void fromFile(const std::string &fn, const std::string &key) {
        std::ifstream in(fn.c_str());
        loadLabels(nw, nh, label, in, key);
        in.close();
    }

    inline void fromString(const std::string &text, const std::string &key) {
        std::istringstream in(text);
        loadLabels(nw, nh, label, in, key);
    }
};

struct Boundary {
    double R;
    double T;
    std::vector<double> unode;
    std::vector<double> dtnode;
    std::vector<bool> where;
    std::string currentVar;

    inline double eval(const std::map<std::string, double> &variables) const {
        double t = variables.at("t");
        t = std::fmod(t, T);
        int iseg = 0;
        double time = 0.0;
        for (int i = 0; i < dtnode.size(); i++) {
            if (time <= t && t < time + dtnode[i]) {
                iseg = i;
                break;
            }
            time += dtnode[i];
        }
        double unodeiseg1;
        if (iseg + 1 == dtnode.size()) unodeiseg1 = unode[0]; else unodeiseg1 = unode[iseg + 1];
        double udrive = (t - time) / dtnode[iseg] * (unodeiseg1 - unode[iseg]) + unode[iseg];

        return udrive - R * variables.at(currentVar);
    }

};

template<typename T>
inline bool lineGetValue(T &val, const std::vector <std::string> &lines,
                         const std::string commandkey, T defval) {

    val = defval;
    for (std::string line : lines) {
        std::istringstream ls(line);
        std::string command;

        ls >> command;
        if (command == commandkey) {
            ls >> val;
            return true;
        }
    }

    return false;

}

inline bool
lineGetSymMatrix(Eigen::MatrixXd &matrix, const std::vector <std::string> &lines, const std::string commandkey,
                 const std::map<std::string, int> &idmap, double defval = 0.0) {

    matrix.resize(idmap.size(), idmap.size());
    matrix.setConstant(defval);
    for (std::string line : lines) {
        std::string command;
        std::istringstream ls(line);
        ls >> command;
        if (command == commandkey) {
            matrix.resize(idmap.size(), idmap.size());
            int no;
            std::string keys1, keys2;
            double value;
            ls >> no;
            for (int k = 0; k < no; k++) {
                ls >> keys1 >> keys2 >> value;

                for (auto key1 : split(keys1, ','))
                    for (auto key2 : split(keys2, ',')) {
                        matrix(idmap.at(key1), idmap.at(key2)) = value;
                        matrix(idmap.at(key2), idmap.at(key1)) = value;
                    }
            }
            return true;
        }
    }
    return false;

}

inline bool lineGetVector(Eigen::VectorXd &vec, const std::vector <std::string> &lines,
                          const std::string commandkey,
                          const std::map<std::string, int> &idmap, double defval = 0.0) {

    vec.resize(idmap.size());
    vec.setConstant(defval);
    for (std::string line : lines) {

        std::string command;

        std::istringstream ls(line);
        ls >> command;
        if (command == commandkey) {

            int no;
            std::string keys1;
            double value;
            ls >> no;
            for (int k = 0; k < no; k++) {
                ls >> keys1 >> value;

                for (auto key1 : split(keys1, ',')) {
                    vec(idmap.at(key1)) = value;
                }
            }
            return true;
        }
    }
    return false;

}

inline bool lineGetIdMap(std::map<std::string, int> &idmap, const std::vector <std::string> &lines,
                         const std::string commandkey) {
    for (std::string line : lines) {
        std::string command;
        std::istringstream ls(line);
        ls >> command;
        if (command == commandkey) {
            int no;
            std::string key;
            ls >> no;
            for (int k = 0; k < no; k++) {
                ls >> key;
                idmap[key] = k;
            }
            return true;
        }
    }
    return false;
}

inline bool lineGetBoundary(std::vector <Boundary> &boundary,
                            const std::vector <std::string> &lines,
                            const std::map<std::string, int> &idmap) {

    boundary.clear();
    for (std::string line : lines) {
        std::istringstream ls(line);
        std::string command;
        ls >> command;
        if (command == "boundary") {

            int no;
            std::string keys1;
            ls >> no;
            for (int k = 0; k < no; k++) {
                Boundary b;
                int numof = 0;
                ls >> keys1 >> b.R >> b.currentVar >> numof;
                b.unode.resize(numof);
                b.dtnode.resize(numof);
                for (int z = 0; z < numof; z++) ls >> b.unode[z] >> b.dtnode[z];
                b.T = 0.0;
                for (int z = 0; z < numof; z++) b.T += b.dtnode[z];
                b.where.resize(idmap.size(), false);
                for (auto key1 : split(keys1, ',')) b.where[idmap.at(key1)] = true;

                boundary.push_back(b);
            }
            return true;
        }
    }
    return false;
}

inline bool lineGetRedox(std::vector < std::tuple < std::vector < bool > , int, int, double, double,
                         double >> &oxidation,
                         std::vector < std::tuple < std::vector < bool > , int, int, double, double,
                         double >> &reduction,
const std::vector <std::string> &lines,
const std::map<std::string, int> &idmap
)
{
oxidation.

clear();

reduction.

clear();

bool oxfound, redfound;
for(
std::string line
: lines) {
std::istringstream ls(line);
std::string command;
if (command == "red") {
int no;
ls >>
no;

reduction.

clear();

for(
int k = 0;
k<no;
k++) {
std::string keys1, key2, key3;
double value1, value2, value3;
ls >> keys1>>key2>>key3 >> value1 >> value2 >>
value3; // where, what, what-to, Ef, dos, tau

std::vector<bool> where(idmap.size(), false);

for(
auto key1
:
split(keys1,
',')) where[idmap.
at(key1)
] = true;

reduction.
push_back(std::tuple<std::vector < bool>,
int,int,double,double,double>(where,idmap.
at(key2), idmap
.
at(key3), value1, value2, value3
));
}
redfound = true;
break;

}
}
for(
std::string line
: lines) {
std::istringstream ls(line);
std::string command;
if (command == "ox") {
int no;
ls >>
no;

for(
int k = 0;
k<no;
k++) {
std::string keys1, key2, key3;
double value1, value2, value3;
ls >> keys1>>key2>>key3 >> value1 >> value2 >>
value3; // where, what, what-to, Ef, dos, tau

std::vector<bool> where(idmap.size(), false);

for(
auto key1
:
split(keys1,
',')) where[idmap.
at(key1)
] = true;
oxidation.
push_back(std::tuple<std::vector < bool>,
int,int,double,double,double>(where,idmap.
at(key2), idmap
.
at(key3), value1, value2, value3
));
}

oxfound = true;
break;

}
}
return
oxfound &&redfound;
}

inline bool lineGetCurrent(std::vector < std::tuple < std::string, std::vector < bool > ,
                           std::vector < bool >> > &current_surfaces,
const std::vector <std::string> &lines,
const std::map<std::string, int> &idmap
)
{
current_surfaces.

clear();

for(
std::string line
: lines) {
std::istringstream ls(line);
std::string command;
ls >>
command;
if (command == "current") {
int no;
ls >>
no;
for(
int k = 0;
k<no;
k++) {
std::string tkeys1;
std::string tkeys2;
std::string name;

ls >> name >> tkeys1 >>
tkeys2;

std::vector<bool> where1(idmap.size(), false);
std::vector<bool> where2(idmap.size(), false);
for(
auto key1
:
split(tkeys1,
',')) where1[idmap.
at(key1)
] = true;
for(
auto key2
:
split(tkeys2,
',')) where2[idmap.
at(key2)
] = true;

current_surfaces.
push_back(std::tuple<std::string, std::vector < bool>, std::vector<bool>>
(name,where1,where2));
}
return true;
}
}
return false;
}
inline bool lineGetVoltage(std::vector < std::tuple < std::string, std::vector < bool >> > &voltage_sites,
const std::vector <std::string> &lines,
const std::map<std::string, int> &idmap
)
{
voltage_sites.

clear();

for(
std::string line
: lines) {
std::istringstream ls(line);
std::string command;
ls >>
command;
if (command == "voltage") {
int no;
ls >>
no;
for(
int k = 0;
k<no;
k++) {
std::string tkeys1;
std::string name;

ls >> name >>
tkeys1;

std::vector<bool> where1(idmap.size(), false);
for(
auto key1
:
split(tkeys1,
',')) where1[idmap.
at(key1)
] = true;

voltage_sites.
push_back(std::tuple<std::string, std::vector < bool>>
(name,where1));
}
return true;
}
}
return false;
}

inline bool lineGetPattern(std::vector <std::tuple<std::vector < std::tuple < int, char, int>>, int, int, double

>> &pattern_transition,
const std::vector <std::string> &lines,
const std::map<std::string, int> &idmap
)
{
pattern_transition.

clear();

for(
std::string line
: lines) {
std::istringstream ls(line);
std::string command;

ls >>
command;
if (command == "pattern") {
int no;
ls >>
no;

pattern_transition.

clear();

for(
int k = 0;
k<no;
k++) {
std::string tkeys1;
std::string tkeys2;
std::string exprs;
double prob;

ls >> exprs >> tkeys1 >> tkeys2 >>
prob;

for(
auto tkey1
:
split(tkeys1,
','))
for(
auto tkey2
:
split(tkeys2,
',')) {

std::vector <std::tuple<int, char, int>> exprli;
for(
auto expr
:
split(exprs,
',')) {

int numof = 0;
std::stringstream ss;
ss << expr[0];
ss >>
numof;
int tt = idmap.at(expr.substr(2));

if( expr[1]=='=' || expr[1]=='<' || expr[1]=='>' ) {
exprli.
push_back( std::tuple<int, char, int>(numof, expr[1], tt)
);
} else throw std::runtime_error("Error in pattern!");
}

pattern_transition.
push_back(std::tuple<std::vector < std::tuple < int, char, int>>
,int,int,double>(exprli,idmap.
at(tkey1), idmap
.
at(tkey2), prob
));

}
}

return true;
}
}

return false;
}

struct Parameters {
    double dt;
    double kT;
    double a;
    double h;
    double omega;
    std::map<std::string, int> lab2type, lab2site;
    std::map<int, std::string> type2lab, site2lab;

    Eigen::MatrixXd ter, tsigma;
    Eigen::VectorXd tsource;
    Eigen::VectorXd tcharge;

    Eigen::MatrixXcd tec;

    std::vector <Boundary> boundary;

    Eigen::VectorXd ttau;
    Eigen::MatrixXd tbarrier;
    Eigen::MatrixXd tgamma1st;
    Eigen::MatrixXd tgamma2nd;

    std::vector<std::tuple<std::vector < bool>, int, int, double, double, double>> oxidation,
    reduction;
    std::vector <std::tuple<std::vector < std::tuple < int, char, int>>,int,int,double>>
    pattern_transition;

    std::vector <std::tuple<std::string, std::vector < bool>, std::vector<bool>>>
    current_surfaces;
    std::vector <std::tuple<std::string, std::vector < bool>>>
    voltage_sites;


    //std::map<std::string,double> variables;


    inline void logStream(std::ostream &os) {
        os << "# dt: " << std::endl << dt << std::endl;
        os << "# a: " << std::endl << a << std::endl;
        os << "# h: " << std::endl << h << std::endl;
        os << "# kT: " << std::endl << kT << std::endl;
        os << "# omega: " << std::endl << omega << std::endl;


        os << "# tcharge: " << std::endl << tcharge.transpose() << std::endl;
        os << "# ttau: " << std::endl << ttau.transpose() << std::endl;
        os << "# tsource: " << std::endl << tsource.transpose() << std::endl;

        os << "# ter: " << std::endl << ter << std::endl;
        os << "# tsigma: " << std::endl << tsigma << std::endl;
        os << "# tec: " << std::endl << tec << std::endl;

        os << "# tgamma1st: " << std::endl << tgamma1st << std::endl;
        os << "# tgamma2nd: " << std::endl << tgamma2nd << std::endl;
        os << "# tbarrier: " << std::endl << tbarrier << std::endl;

        os << "# lab2type: " << std::endl;
        for (auto it = lab2type.begin(); it != lab2type.end(); it++) os << it->first << ':' << it->second << '\t';
        os << std::endl;

        os << "# lab2site: " << std::endl;
        for (auto it = lab2site.begin(); it != lab2site.end(); it++) os << it->first << ':' << it->second << '\t';
        os << std::endl;

        os << "# type2lab: " << std::endl;
        for (auto it = type2lab.begin(); it != type2lab.end(); it++) os << it->first << ':' << it->second << '\t';
        os << std::endl;

        os << "# site2lab: " << std::endl;
        for (auto it = site2lab.begin(); it != site2lab.end(); it++) os << it->first << ':' << it->second << '\t';
        os << std::endl;

        os << "# oxidation: " << std::endl;
        for (auto e : oxidation) {
            for (bool b : std::get<0>(e)) os << b << ' ';
            os << '\t';
            os << std::get<1>(e) << '\t';
            os << std::get<2>(e) << '\t';
            os << std::get<3>(e) << '\t';
            os << std::get<4>(e) << '\t';
            os << std::get<5>(e) << '\t';
            os << std::endl;
        }

        os << "# reduction: " << std::endl;
        for (auto e : reduction) {
            for (bool b : std::get<0>(e)) os << b << ' ';
            os << '\t';
            os << std::get<1>(e) << '\t';
            os << std::get<2>(e) << '\t';
            os << std::get<3>(e) << '\t';
            os << std::get<4>(e) << '\t';
            os << std::get<5>(e) << '\t';
            os << std::endl;
        }

        os << "# pattern_transition: " << std::endl;
        for (auto e : pattern_transition) {
            for (auto ici : std::get<0>(e))
                os << std::get<0>(ici) << ' ' << std::get<1>(ici) << ' ' << std::get<2>(ici) << ' ';
            os << '\t';
            os << std::get<1>(e) << '\t';
            os << std::get<2>(e) << '\t';
            os << std::get<3>(e) << '\t';
            os << std::endl;
        }

        os << "# current_surfaces: " << std::endl;
        for (auto e : current_surfaces) {
            os << std::get<0>(e) << '\t';
            for (bool b : std::get<1>(e)) os << b << ' ';
            os << '\t';
            for (bool b : std::get<2>(e)) os << b << ' ';
            os << '\t';

            os << std::endl;
        }
        os << "# voltage_sites: " << std::endl;
        for (auto e : voltage_sites) {
            os << std::get<0>(e) << '\t';
            for (bool b : std::get<1>(e)) os << b << ' ';
            os << '\t';

            os << std::endl;
        }

        os << "# boundary: " << std::endl;
        for (auto b : boundary) {

            for (bool bo : b.where) os << bo << ' ';
            os << '\t';
            os << b.T << '\t';
            os << b.R << '\t';
            os << b.currentVar << '\t';
            for (auto uu : b.unode) os << uu << ' ';
            os << '\t';
            for (auto dt : b.dtnode) os << dt << ' ';
            os << '\t';

            os << std::endl;
        }


    }


    inline void fromStream(std::istream &is) {

        std::vector <std::string> lines;
        std::string line;
        while (std::getline(is, line)) lines.push_back(line);


        lineGetIdMap(lab2type, lines, "type");
        lineGetIdMap(lab2site, lines, "site");

        lineGetValue(dt, lines, "dt", 1.0);
        lineGetValue(kT, lines, "kT", 1.0);
        lineGetValue(a, lines, "a", 1.0);
        lineGetValue(h, lines, "h", 1.0);
        lineGetValue(omega, lines, "omega", 0.001);

        lineGetVector(ttau, lines, "tau", lab2type, 1.0);
        lineGetVector(tsource, lines, "source", lab2type, 0.0);
        lineGetVector(tcharge, lines, "charge", lab2type, 0.0);
        lineGetSymMatrix(ter, lines, "er", lab2type, 1.0);
        lineGetSymMatrix(tsigma, lines, "sigma", lab2type, 0.0);
        lineGetSymMatrix(tbarrier, lines, "barrier", lab2type, 0.0);
        lineGetSymMatrix(tgamma1st, lines, "gamma1st", lab2type, 0.0);
        lineGetSymMatrix(tgamma2nd, lines, "gamma2nd", lab2type, 0.0);

        lineGetBoundary(boundary, lines, lab2site);

        lineGetRedox(oxidation, reduction, lines, lab2type);
        lineGetPattern(pattern_transition, lines, lab2type);

        lineGetCurrent(current_surfaces, lines, lab2site);
        lineGetVoltage(voltage_sites, lines, lab2site);

        for (auto it = lab2type.begin(); it != lab2type.end(); it++) type2lab[it->second] = it->first;
        for (auto it = lab2site.begin(); it != lab2site.end(); it++) site2lab[it->second] = it->first;

        tec.resize(lab2type.size(), lab2type.size());
        tec = std::complex<double>(1.0, 0.0) * ter - 4.0 * 3.1415 / omega * std::complex<double>(0.0, 1.0) * tsigma;
    }

    inline void fromFile(const std::string &fn) {
        std::ifstream in(fn.c_str());
        fromStream(in);
        in.close();
    }

    inline void fromString(const std::string &text) {
        std::istringstream in(text);
        fromStream(in);
    }

    inline void gettypes(Eigen::VectorXi &types, const std::vector <std::string> &labels) const {
        types.resize(labels.size());
        for (int k = 0; k < labels.size(); k++) types(k) = lab2type.at(labels[k]);
    }

    inline void getsites(Eigen::VectorXi &sites, const std::vector <std::string> &labels) const {
        sites.resize(labels.size());
        for (int k = 0; k < labels.size(); k++) sites(k) = lab2site.at(labels[k]);
    }

    inline void gettypelabel(std::vector <std::string> &labels, const Eigen::VectorXi &types) const {
        labels.resize(types.rows());
        for (int k = 0; k < labels.size(); k++) labels[k] = type2lab.at(types(k));
    }

    inline void getsitelabel(std::vector <std::string> &labels, const Eigen::VectorXi &sites) const {
        labels.resize(sites.rows());
        for (int k = 0; k < labels.size(); k++) labels[k] = site2lab.at(sites(k));
    }

    inline void getec(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> &ec,
                      const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb, const Eigen::VectorXi &types) const {
        ec.resize(nb.rows(), Eigen::NoChange);
        for (int k = 0; k < nb.rows(); k++) {
            for (int i = 0; i < 6; i++)
                ec(k, i) = tec(types(k), types(nb(k, i)));
        }
    }

    inline void getsigma(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> &sigma,
                         const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb, const Eigen::VectorXi &types) const {
        sigma.resize(nb.rows(), Eigen::NoChange);
        for (int k = 0; k < nb.rows(); k++) {
            for (int i = 0; i < 6; i++)
                sigma(k, i) = tsigma(types(k), types(nb(k, i)));
        }
    }

    inline void getsource(Eigen::Matrix<double, Eigen::Dynamic, 1> &source, const Eigen::VectorXi &types) const {
        source.resize(types.size());
        for (int k = 0; k < types.size(); k++) source(k) = tsource(types(k));

    }

    inline void getboundary(Eigen::VectorXi &B, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &u0,
                            const Eigen::VectorXi &sites, const std::map<std::string, double> &variables) const {
        Eigen::VectorXi sboundary(lab2site.size());
        sboundary.setZero();
        Eigen::VectorXd su0(lab2site.size());
        su0.setZero();
        for (int s = 0; s < lab2site.size(); s++)
            for (int i = 0; i < boundary.size(); i++)
                if (boundary[i].where[s]) {
                    sboundary(s) = 1;
                    su0(s) = boundary[i].eval(variables);
                }
        B.resize(sites.rows());
        u0.resize(sites.rows());
        for (int k = 0; k < sites.size(); k++) {
            B(k) = sboundary(sites(k));
            u0(k) = su0(sites(k));
        }
    }

    inline void geteb(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> &e,
                      Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &b,
                      const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                      const Eigen::VectorXi &types, const Eigen::VectorXi &sites,
                      const std::map<std::string, double> &variables) const {

        Eigen::VectorXi sboundary(lab2site.size());
        sboundary.setZero();
        Eigen::VectorXd su0(lab2site.size());
        su0.setZero();
        for (int s = 0; s < lab2site.size(); s++)
            for (int i = 0; i < boundary.size(); i++)
                if (boundary[i].where[s]) {
                    sboundary(s) = 1;
                    su0(s) = boundary[i].eval(variables);
                }

        e.resize(nb.rows(), 6);
        b.resize(nb.rows());

        const double sourceFactor = 4.0 * 3.1415 * std::sqrt(3.0) / h;

        for (int k = 0; k < e.rows(); k++) {
            int sk = sites(k);

            if (sboundary(sk) == 1) {
                e.row(k).setZero();
                b(k) = su0(sk);
            } else {
                int tk = types(k);

                for (int i = 0; i < 6; i++) e(k, i) = tec(tk, types(nb(k, i)));
                std::complex<double> isum = 1.0 / e.row(k).sum();

                e.row(k) *= isum;
                b(k) = tsource(tk) * sourceFactor * isum;
            }

        }

    }

    inline void getIEq(Eigen::Matrix<double, Eigen::Dynamic, 6> &I,
                       Eigen::Matrix<double, Eigen::Dynamic, 6> &E,
                       Eigen::Matrix<double, Eigen::Dynamic, 1> &q,
                       const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &u,
                       const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                       const Eigen::VectorXi &types) const {

        I.resize(nb.rows(), 6);
        E.resize(nb.rows(), 6);
        q.resize(nb.rows());

        const double sigmaFactor = h / std::sqrt(3.0);

        for (int k = 0; k < nb.rows(); k++) {
            q(k) = 0.0;
            for (int i = 0; i < 6; i++) {
                double ski = tsigma(types(k), types(nb(k, i)));
                E(k, i) = -std::real(u(nb(k, i)) - u(k)) / a;
                I(k, i) = -std::real(u(nb(k, i)) - u(k)) * sigmaFactor * ski;
                q(k) += std::imag(u(nb(k, i)) - u(k)) * sigmaFactor * ski / omega;

            }
        }

    }

    inline void getVarCurrent(std::map<std::string, double> &variables,
                              const Eigen::Matrix<double, Eigen::Dynamic, 6> &I,
                              const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                              const Eigen::VectorXi &sites) const {

        for (auto e : current_surfaces) {

            double current = 0.0;
            auto where1 = std::get<1>(e);
            auto where2 = std::get<2>(e);
            for (int k = 0; k < I.rows(); k++)
                if (where1[sites(k)]) {
                    for (int i = 0; i < 6; i++)
                        if (where2[sites(nb(k, i))]) {
                            current += I(k, i);
                        }
                }
            variables[std::get<0>(e)] = current;
        }


    }

    inline void getVarVoltage(std::map<std::string, double> &variables,
                              const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
                              const Eigen::VectorXi &sites) const {

        for (auto e : voltage_sites) {
            int nosites = 0;
            double voltage = 0.0;
            auto where1 = std::get<1>(e);
            for (int k = 0; k < u.rows(); k++)
                if (where1[sites(k)]) {
                    voltage += u(k);
                    nosites++;
                }

            variables[std::get<0>(e)] = nosites > 0 ? voltage / nosites : 0.0;
        }


    }

    inline void gettypenb(Eigen::Matrix<int, Eigen::Dynamic, 13> &typenb,
                          const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                          const Eigen::VectorXi &types) const {
        typenb = nb;
        for (int r = 0; r < typenb.rows(); r++)
            for (int c = 0; c < typenb.cols(); c++)
                typenb(r, c) = types(nb(r, c));
    }

    inline double chempot(int k, int t, const Eigen::Matrix<int, Eigen::Dynamic, 13> &typenb) const {
        double mu = 0.0;
        for (int i = 0; i < 6; i++) mu += tgamma1st(t, typenb(k, i));
        for (int i = 6; i < 12; i++) mu += tgamma2nd(t, typenb(k, i));
        return mu;
    }

    inline double elchempot(int k, int t, const Eigen::Matrix<int, Eigen::Dynamic, 13> &typenb,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &u) const {
        return chempot(k, t, typenb) + u(k) * tcharge(t);
    }

    inline void diffusion(Eigen::VectorXi &types,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
                          const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb_array,
                          std::default_random_engine &generator) const {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        std::vector <std::pair<int, int>> transitions;
        transitions.reserve(4 * types.rows());
        for (int a = 0; a < types.rows(); a++)
            for (int nbi = 0; nbi < 6; nbi++) {
                int b = nb_array(a, nbi);
                if (a < b)
                    if (types(a) != types(b)) {
                        transitions.push_back(std::pair<int, int>(a, b));
                    }
            }

        std::shuffle(transitions.begin(), transitions.end(), generator);

        int nofswap = 0;
        double maxee = 0.0;
        double minee = 0.0;
        for (int i = 0; i < transitions.size(); i++) {
            int a = transitions[i].first;
            int b = transitions[i].second;
            int ta = types[a];
            int tb = types[b];

            if (ta == tb) continue;

            double taua = ttau(ta);
            double taub = ttau(tb);

            //potential charge
            double eabu = u(a) * tcharge(tb) + u(b) * tcharge(ta) - (u(a) * tcharge(ta) + u(b) * tcharge(tb));
            //  pair interaction
            double eab = 0.0;
            for (int ni = 0; ni < 6; ni++) {
                int na = nb_array(a, ni);
                eab += tgamma1st(tb, types[na]) - tgamma1st(ta, types[na]);
                int nb = nb_array(b, ni);
                eab += tgamma1st(ta, types[nb]) - tgamma1st(tb, types[nb]);
            }
            eab += 2.0 * tgamma1st(ta, tb) - tgamma1st(ta, ta) - tgamma1st(tb, tb);

            //  2nd nb pair interaction
            double eab2nd = 0.0;
            for (int ni = 6; ni < 12; ni++) {
                int na = nb_array(a, ni);
                eab2nd += tgamma2nd(tb, types[na]) - tgamma2nd(ta, types[na]);
                int nb = nb_array(b, ni);
                eab2nd += tgamma2nd(ta, types[nb]) - tgamma2nd(tb, types[nb]);
            }

            // barrier
            double eabb = tbarrier(ta, tb);

            double ee = eabu + eab + eab2nd + eabb;

            double w = std::min(1.0, dt / (taua + taub) * std::exp(-ee / kT));

            if (distribution(generator) < w) {
                types[a] = tb;
                types[b] = ta;
                nofswap++;
            }

            maxee = std::max(maxee, ee);
            minee = std::min(minee, ee);

        }
    }

    inline void redox(Eigen::VectorXi &types,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &qs,
                      const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb_array,
                      std::default_random_engine &generator) const {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);


        Eigen::Matrix<int, Eigen::Dynamic, 13> typenb;
        gettypenb(typenb, nb_array, types);

        std::vector <std::tuple<int, int>> red_transition;
        for (int r = 0; r < reduction.size(); r++) {
            const std::vector<bool> &where = std::get<0>(reduction[r]);
            int from = std::get<1>(reduction[r]);

            for (int a = 0; a < types.rows(); a++)
                if (types[a] == from) {
                    bool surf = false;
                    for (int nbi = 0; nbi < 6; nbi++)
                        if (where[typenb(a, nbi)]) {
                            surf = true;
                            break;
                        }
                    if (surf) red_transition.push_back(std::tuple<int, int>(a, r));
                }

        }
        std::vector <std::tuple<int, int>> ox_transition;
        for (int r = 0; r < oxidation.size(); r++) {
            const std::vector<bool> &where = std::get<0>(oxidation[r]);
            int from = std::get<1>(oxidation[r]);
            for (int a = 0; a < types.rows(); a++)
                if (types[a] == from) {
                    bool surf = false;
                    for (int nbi = 0; nbi < 6; nbi++)
                        if (where[typenb(a, nbi)]) {
                            surf = true;
                            break;
                        }
                    if (surf) ox_transition.push_back(std::tuple<int, int>(a, r));
                }
        }

        std::vector <std::tuple<int, int>> redoxpairs;
        for (int kr = 0; kr < red_transition.size(); kr++)
            for (int ko = 0; ko < ox_transition.size(); ko++) redoxpairs.push_back(std::tuple<int, int>(kr, ko));

        std::shuffle(redoxpairs.begin(), redoxpairs.end(), generator);

        int nofswap = 0;
        double maxee = 0.0;
        double minee = 0.0;

        for (int ro = 0; ro < redoxpairs.size(); ro++) {

            int kr = std::get<0>(redoxpairs[ro]);
            int ko = std::get<1>(redoxpairs[ro]);

            int ari = std::get<0>(red_transition[kr]);
            int rr = std::get<1>(red_transition[kr]);
            int tarifr = std::get<1>(reduction[rr]);
            int tarito = std::get<2>(reduction[rr]);
            double efred = std::get<3>(reduction[rr]);
            double dosred = std::get<4>(reduction[rr]);
            double taured = std::get<5>(reduction[rr]);

            int aoi = std::get<0>(ox_transition[ko]);
            int oo = std::get<1>(ox_transition[ko]);
            int taoifr = std::get<1>(oxidation[oo]);
            int taoito = std::get<2>(oxidation[oo]);
            double efox = std::get<3>(oxidation[oo]);
            double dosox = std::get<4>(oxidation[oo]);
            double tauox = std::get<5>(oxidation[oo]);

            if (tarifr != types[ari] || taoifr != types[aoi]) continue;

            double muox = efox - u(aoi) - qs(aoi) / dosox;
            double dEox = elchempot(aoi, taoito, typenb, u) - elchempot(aoi, taoifr, typenb, u) + muox;

            double mured = efred - u(ari) - qs(ari) / dosred;
            double dEred = elchempot(ari, tarito, typenb, u) - elchempot(ari, tarifr, typenb, u) - mured;

            double ee = dEox + dEred;
            double w = std::min(1.0, dt / (tauox + taured) * std::exp(-ee / kT));

            if (distribution(generator) < w) {
                types[ari] = tarito;
                types[aoi] = taoito;
                nofswap++;
            }

            maxee = std::max(maxee, ee);
            minee = std::min(minee, ee);
        }
    }

    inline void pattern(Eigen::VectorXi &types,
                        const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb_array,
                        std::default_random_engine &generator) const {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        auto temptypes = types;
        for (int a = 0; a < types.rows(); a++) {
            for (int l = 0; l < pattern_transition.size(); l++) {

                auto cond = std::get<0>(pattern_transition[l]);
                int tfr = std::get<1>(pattern_transition[l]);
                int tto = std::get<2>(pattern_transition[l]);
                double prob = std::get<3>(pattern_transition[l]);

                if (tfr == types[a]) {

                    bool condsatistfy = true;
                    for (int c = 0; c < cond.size(); c++) {
                        int numof = std::get<0>(cond[c]);
                        char cha = std::get<1>(cond[c]);
                        int typn = std::get<2>(cond[c]);

                        int counter = 0;
                        for (int ni = 0; ni < 6; ni++) if (types[nb_array(a, ni)] == typn) counter++;

                        if (cha == '=') {
                            condsatistfy = condsatistfy && (numof == counter);
                        } else if (cha == '<') {
                            condsatistfy = condsatistfy && (numof < counter);
                        } else if (cha == '>') {
                            condsatistfy = condsatistfy && (numof > counter);
                        }
                    }
                    if (cond.size() < 0) condsatistfy = false;

                    if (condsatistfy && distribution(generator) <= prob) {
                        temptypes[a] = tto;
                    }
                }


            }
        }

        types = temptypes;

    }


};

struct Control {

    std::default_random_engine generator;

    Eigen::Matrix<int, Eigen::Dynamic, 13> nb;
    Eigen::Matrix<double, Eigen::Dynamic, 2> co;
    Eigen::Matrix<double, 6, 2> an;

    Eigen::VectorXi types, sites;
    int nw, nh;

    Eigen::VectorXi B;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> u0;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> u;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> e;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> b;

    Eigen::Matrix<double, Eigen::Dynamic, 6> I;
    Eigen::Matrix<double, Eigen::Dynamic, 6> E;
    Eigen::Matrix<double, Eigen::Dynamic, 1> q;

    bool flagPoisson;
    bool flagDiffusion;
    bool flagRedox;
    bool flagPattern;

    std::map<std::string, double> variables;

    inline void initVariables() {

        variables["_"] = 0.0;
        variables["t"] = 0.0;
        variables["$1"] = 0.0;
        variables["$2"] = 0.0;
        variables["$3"] = 0.0;
        variables["$4"] = 0.0;
        variables["$5"] = 0.0;
        variables["$6"] = 0.0;
        variables["$7"] = 0.0;
        variables["$8"] = 0.0;
        variables["$9"] = 0.0;
        variables["$0"] = 0.0;
    }

    inline void setup(const Label &state, const Label &site, const Parameters &parameters) {
        assert(state.nw == site.nw && state.nh == site.nh);
        nw = state.nw;
        nh = state.nh;

        generateGrid(nb, co, an, state.nw, state.nh);

        parameters.gettypes(types, state.label);
        parameters.getsites(sites, site.label);

        parameters.getboundary(B, u0, sites, variables);

        parameters.geteb(e, b, nb, types, sites, variables);

        sparseMatrixSolver(u, nb, e, b);
        parameters.getIEq(I, E, q, u, nb, types);

        parameters.getVarCurrent(variables, I, nb, sites);

        parameters.getVarVoltage(variables, u.real(), sites);
    }

    void getCurrentState(Label &getstate, const Label &state, const Parameters &parameters) {
        getstate = state;
        parameters.gettypelabel(getstate.label, types);
    }

    inline void step(const Parameters &parameters) {
        if (flagDiffusion) parameters.diffusion(types, u.real(), nb, generator);
        if (flagRedox) parameters.redox(types, u.real(), q, nb, generator);
        if (flagPattern) parameters.pattern(types, nb, generator);
        if (flagPoisson) {
            parameters.geteb(e, b, nb, types, sites, variables);
            sparseMatrixSolver(u, nb, e, b);
            parameters.getIEq(I, E, q, u, nb, types);
            parameters.getVarCurrent(variables, I, nb, sites);
            parameters.getVarVoltage(variables, u.real(), sites);
        }

        variables["t"] += parameters.dt;
    }

    inline void logStream(std::ostream &os) {

        os << variables["t"] << '\t';
        for (auto v : variables) if (v.first != "_" && v.first != "t") os << v.second << '\t';
        os << std::endl;
        os << nw << ' ' << nh;
        os << std::endl;
        os << u.real().transpose() << std::endl;
        os << std::endl;
        os << q.transpose() << std::endl;
        os << std::endl;
        os << types.transpose() << std::endl;
    }


    int run;
    int every;

    inline void fromStream(std::istream &is) {

        std::vector <std::string> lines;
        std::string line;
        while (std::getline(is, line)) lines.push_back(line);

        int seed;
        lineGetValue(seed, lines, "seed", 0);
        generator.seed(seed);

        lineGetValue(run, lines, "run", 1);
        lineGetValue(every, lines, "every", 1);


        lineGetValue(flagPoisson, lines, "flagPoisson", true);
        lineGetValue(flagDiffusion, lines, "flagDiffusion", false);
        lineGetValue(flagRedox, lines, "flagRedox", false);
        lineGetValue(flagPattern, lines, "flagPattern", false);


        double time;
        if (lineGetValue(time, lines, "time", 0.0)) {
            variables["t"] = time;
        }
    }

    inline void fromFile(const std::string &fn) {
        std::ifstream in(fn.c_str());
        fromStream(in);
        in.close();
    }

    inline void fromString(const std::string &text) {
        std::istringstream in(text);
        fromStream(in);
    }
};

