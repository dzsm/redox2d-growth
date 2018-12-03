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

#ifndef GROWTH_BASE_H
#define GROWTH_BASE_H

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


template<typename T>
void sparseMatrixSolver(Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
                        const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                        const Eigen::Matrix<T, Eigen::Dynamic, 6> &e,
                        const Eigen::Matrix<T, Eigen::Dynamic, 1> &b) {
    if (b.rows() == 0 || e.rows() == 0 || nb.rows() == 0) return;

    if (u.rows() != b.rows()) {
        u.resize(b.rows());
        u.setZero();
    }

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

struct SparseMatrixSolver {

    Eigen::BiCGSTAB <Eigen::SparseMatrix<std::complex < double>>,Eigen::IncompleteLUT <std::complex<double>>>
    solver;

    void solveFirst(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &u,
                    const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
                    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> &e,
                    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &b) {
        if (b.rows() == 0 || e.rows() == 0 || nb.rows() == 0) {
            return;
        }

        if (u.rows() != b.rows()) {
            u.resize(b.rows());
            u.setZero();
        }

        Eigen::SparseMatrix <std::complex<double>> poisson;        // default is column major

        poisson.resize(nb.rows(), nb.rows());         // default is column major
        poisson.reserve(Eigen::VectorXi::Constant(nb.rows(), 7));

        for (int k = 0; k < nb.rows(); k++) {
            for (int i = 0; i < 6; i++)
                poisson.insert(k, nb(k, i)) = e(k, i);
            poisson.insert(k, k) = -1.0;
        }
        poisson.makeCompressed();
        solver.analyzePattern(poisson);
        solver.preconditioner().compute(poisson);
        solver.factorize(poisson);
        u = solver.solveWithGuess(-b, u);
    }

    void solve(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &u,
               const Eigen::Matrix<int, Eigen::Dynamic, 13> &nb,
               const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> &e,
               const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &b) {
        assert(b.rows() == e.rows() && e.rows() == nb.rows() && nb.rows() == u.rows());
        if (b.rows() == 0 || e.rows() == 0 || nb.rows() == 0) {
            return;
        }

        Eigen::SparseMatrix <std::complex<double>> poisson;        // default is column major

        poisson.resize(nb.rows(), nb.rows());         // default is column major
        poisson.reserve(Eigen::VectorXi::Constant(nb.rows(), 7));

        for (int k = 0; k < nb.rows(); k++) {
            for (int i = 0; i < 6; i++)
                poisson.insert(k, nb(k, i)) = e(k, i);
            poisson.insert(k, k) = -1.0;
        }
        poisson.makeCompressed();

        /*for (int k=0; k<poisson.outerSize(); ++k)
            for (SparseMatrix<std::complex<double>>::InnerIterator it(poisson,k); it; ++it) {
                it.value();
                it.row();   // row index
                it.col();   // col index (here it is equal to k)
                it.index(); // inner index, here it is equal to it.row()
            }
        */
        solver.preconditioner().compute(poisson);

        solver.factorize(poisson);
        u = solver.solveWithGuess(-b, u);
    }
};

struct MultiGridPoissonSolver {

    typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> TVec;
    typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> TSMat;
    typedef Eigen::Matrix<int, Eigen::Dynamic, 6> TIMat;

    double residual_number;

    int MG_PRESMOOTH;
    int MG_POSTSMOOTH;
    int MG_FULL_NCYCLES_FINE;
    double MG_FINAL_TOLERANCE;
    bool doBICG;
    double alpha;

    MultiGridPoissonSolver() {
        alpha = 1.6;
        MG_PRESMOOTH = 3;
        MG_POSTSMOOTH = 3;
        MG_FULL_NCYCLES_FINE = 10;
        MG_FINAL_TOLERANCE = 1e-5;
        doBICG = true;
    }

    void gen_nb6(TIMat &nb6, int nw, int nh) {
        nb6.resize(nw * nh, 6);
        for (int j = 0; j < nh; j++)
            for (int i = 0; i < nw; i++) {
                int k = j * nw + i;
                if (j % 2 == 0) {
                    nb6(k, 0) = index(i + 1, j, nw, nh);
                    nb6(k, 1) = index(i, j + 1, nw, nh);
                    nb6(k, 2) = index(i - 1, j + 1, nw, nh);
                    nb6(k, 3) = index(i - 1, j, nw, nh);
                    nb6(k, 4) = index(i - 1, j - 1, nw, nh);
                    nb6(k, 5) = index(i, j - 1, nw, nh);
                } else {
                    nb6(k, 0) = index(i + 1, j, nw, nh);
                    nb6(k, 1) = index(i + 1, j + 1, nw, nh);
                    nb6(k, 2) = index(i, j + 1, nw, nh);
                    nb6(k, 3) = index(i - 1, j, nw, nh);
                    nb6(k, 4) = index(i, j - 1, nw, nh);
                    nb6(k, 5) = index(i + 1, j - 1, nw, nh);
                }
            }

    }

    int get_depth(int nw, int nh) {
        if (nw == 0 || nh == 0) return 0;
        const int maxLevelPower = 100;
        int div = 1;
        for (int d = 1; d < maxLevelPower; d++) {
            if (nw % 2 == 0 && nh % 2 == 0) {
                nw = nw / 2;
                nh = nh / 2;
                div++;
            } else {
                break;
            }
        }
        return div;
    }

    std::vector <TVec> xl;
    std::vector <TSMat> Al;
    std::vector <TVec> bl;
    std::vector <TIMat> Il;
    std::vector<int> nwl;
    std::vector<int> nhl;

    void generateLadder(const TVec &x, const TSMat &A, const TVec &b, int nw, int nh) {
        xl.clear();
        Al.clear();
        bl.clear();
        Il.clear();
        nwl.clear();
        nhl.clear();

        int depth = get_depth(nw, nh);
        xl.resize(depth);
        Al.resize(depth);
        bl.resize(depth);
        Il.resize(depth);
        nwl.resize(depth);
        nhl.resize(depth);
        xl[depth - 1] = x;
        Al[depth - 1] = A;
        bl[depth - 1] = b;
        gen_nb6(Il[depth - 1], nw, nh);
        nwl[depth - 1] = nw;
        nhl[depth - 1] = nh;

        for (int l = depth - 1; l > 0; l--) {
            restrict(bl[l - 1], Il[l - 1], bl[l], Il[l], nwl[l], nhl[l]);
            restrict(Al[l - 1], Il[l - 1], Al[l], Il[l], nwl[l], nhl[l]);
            nwl[l - 1] = nwl[l] / 2;
            nhl[l - 1] = nhl[l] / 2;
            gen_nb6(Il[l - 1], nwl[l - 1], nhl[l - 1]);
        }

    }

    void updateLadder(const TSMat &A, const TVec &b) {

        int depth = Al.size();
        Al[depth - 1] = A;
        bl[depth - 1] = b;
        for (int l = depth - 1; l > 0; l--) {
            restrict(bl[l - 1], Il[l - 1], bl[l], Il[l], nwl[l], nhl[l]);
            restrict(Al[l - 1], Il[l - 1], Al[l], Il[l], nwl[l], nhl[l]);
        }

    }


    void solveFirst(TVec &x,
                    const TSMat &A,
                    const TVec &b, int nw, int nh) {
        assert(x.rows() == A.rows() && A.rows() == b.rows() && nw * nh == b.rows());
        assert(get_depth(nw, nh) > 1);
        x.setZero();
        generateLadder(x, A, b, nw, nh);
        FMG();
        x = xl[xl.size() - 1];
    }

    void solve(TVec &x,
               const TSMat &A,
               const TVec &b) {
        updateLadder(A, b);
        FMG();
        x = xl[xl.size() - 1];
    }

    void FMG() {
        assert(xl.size() == Al.size() && Al.size() == bl.size() && bl.size() == Il.size() && Il.size() > 0);

        exact(xl[0], Al[0], bl[0], Il[0]);
        for (int l = 1; l < Al.size(); l++) {
            prolong(xl[l], Il[l], xl[l - 1], Il[l - 1], nwl[l - 1], nhl[l - 1]);
            MGV(xl[l], bl[l], l);
        }
        //for(int iter = 0; iter<MG_FULL_NCYCLES_FINE; iter++)  MGV(xl[xl.size()-1],bl[xl.size()-1],xl.size()-1);

        residual_number = 0.0;
        for (int iter = 0; iter < MG_FULL_NCYCLES_FINE; iter++) {
            TVec r;
            residual(r, xl[xl.size() - 1], Al[xl.size() - 1], bl[xl.size() - 1], Il[xl.size() - 1], nwl[xl.size() - 1],
                     nhl[xl.size() - 1]);
            residual_number = r.norm();
            std::cout << '.';
            if (residual_number < MG_FINAL_TOLERANCE) break;
            MGV(xl[xl.size() - 1], bl[xl.size() - 1], xl.size() - 1);
        }
        std::cout << std::endl;
        std::cout << "Poisson residual: " << residual_number << std::endl;

        if (residual_number > MG_FINAL_TOLERANCE && doBICG) {
            TVec r;
            residual(r, xl[xl.size() - 1], Al[xl.size() - 1], bl[xl.size() - 1], Il[xl.size() - 1], nwl[xl.size() - 1],
                     nhl[xl.size() - 1]);
            residual_number = r.norm();
            exact(xl[xl.size() - 1], Al[xl.size() - 1], bl[xl.size() - 1], Il[xl.size() - 1]);
            std::cout << "Poisson residual: " << residual_number << std::endl;
            std::cout << "Warning: exact poisson had to be used!" << std::endl;

        }
    }

    void MGV(TVec &x, const TVec &b, int l) {

        if (l == 0) {
            exact(x, Al[0], b, Il[0]);
        } else {
            for (int iter = 0; iter < MG_PRESMOOTH; iter++) smooth(x, Al[l], b, Il[l], nwl[l], nhl[l]);
            TVec r, rc;
            residual(r, x, Al[l], b, Il[l], nwl[l], nhl[l]);
            restrict(rc, Il[l - 1], r, Il[l], nwl[l],
                     nhl[l]);     //  r = A.x - x - b ;      A.v - v - r = 0 --> A.(x-v) - (x-v) - b = 0 ;
            TVec vc(rc.rows());
            vc.setZero();
            MGV(vc, rc, l - 1);
            TVec v;
            prolong(v, Il[l], vc, Il[l - 1], nwl[l - 1], nhl[l - 1]);
            x -= v;
            for (int iter = 0; iter < MG_POSTSMOOTH; iter++) smooth(x, Al[l], b, Il[l], nwl[l], nhl[l]);
        }

    }

    // L.x = b
    // A.x - x = b
    // L = A-I
    // r = L.x - b
    // r = A.x - x - b;
    void residual(TVec &r, const TVec &x, const TSMat &A, const TVec &b, const TIMat &I, int nw, int nh) {
        r.resize(x.rows());
        for (int k = 0; k < A.rows(); k++) {
            std::complex<double> sum = -b(k) - x(k);
            for (int c = 0; c < 6; c++) sum += A(k, c) * x(I(k, c));
            r(k) = sum;
        }
    }

    // x' = x + alpha*r
    void smooth(TVec &x, const TSMat &A, const TVec &b, const TIMat &I, int nw, int nh) {
        for (int j = 0; j < nh; j++)
            for (int i = j % 2 + 0; i < nw; i += 3) {
                int k = index(i, j, nw, nh);
                std::complex<double> sum = -b(k) - x(k);
                for (int l = 0; l < 6; l++) sum += x(I(k, l)) * A(k, l);
                x(k) += alpha * sum;
            }
        for (int j = 0; j < nh; j++)
            for (int i = j % 2 + 1; i < nw; i += 3) {
                int k = index(i, j, nw, nh);
                std::complex<double> sum = -b(k) - x(k);
                for (int l = 0; l < 6; l++) sum += x(I(k, l)) * A(k, l);
                x(k) += alpha * sum;
            }
        for (int j = 0; j < nh; j++)
            for (int i = (j % 2 + 2) % 3; i < nw; i += 3) {
                int k = index(i, j, nw, nh);
                std::complex<double> sum = -b(k) - x(k);
                for (int l = 0; l < 6; l++) sum += x(I(k, l)) * A(k, l);
                x(k) += alpha * sum;
            }
    }

    void prolong(TVec &xp, const TIMat &Ip, const TVec &xc, const TIMat &Ic, int nw, int nh) {
        xp.resize(nw * nh * 4);

        int nwc = nw;
        int nhc = nh;
        int nwp = 2 * nw;
        int nhp = 2 * nh;

        for (int j = 0; j < nhc; j++)
            for (int i = 0; i < nwc; i++) {
                int ck = index(i, j, nwc, nhc);
                int pk = index(2 * i + j % 2, 2 * j, nwp, nhp);

                xp(pk) = xc(ck);
                xp(Ip(pk, 0)) = 0.5 * (xc(ck) + xc(Ic(ck, 0)));
                xp(Ip(pk, 1)) = 0.5 * (xc(ck) + xc(Ic(ck, 1)));
                xp(Ip(Ip(pk, 0), 1)) = 0.5 * (xc(Ic(ck, 0)) + xc(Ic(ck, 1)));
            }
    }

    void restrict(TVec &xc, const TIMat &Ic, const TVec &x, const TIMat &I, int nw, int nh) {
        xc.resize(nw * nh / 4);
        for (int j = 0; j < nh / 2; j++)
            for (int i = 0; i < nw / 2; i++) {
                int ck = index(i, j, nw / 2, nh / 2);
                int k = index(2 * i + j % 2, 2 * j, nw, nh);
                xc(ck) = x(k);
                //std::cout << ck << ':' << k << ' ';
            }
        // std::cout << "----------" << std::endl;
    }

    void restrict(TSMat &Ac, const TIMat &Ic, const TSMat &A, const TIMat &I, int nw, int nh) {
        Ac.resize(nw * nh / 4, 6);
        for (int j = 0; j < nh / 2; j++)
            for (int i = 0; i < nw / 2; i++) {
                int ck = index(i, j, nw / 2, nh / 2);
                int k = index(2 * i + j % 2, 2 * j, nw, nh);
                for (int l = 0; l < 6; l++) Ac(ck, l) = 0.5 * (A(k, l) + A(I(k, l), l));
            }
    }

    void exact(TVec &x, const TSMat &A, const TVec &b, const TIMat &I) {

        Eigen::SparseMatrix <std::complex<double>> poisson(I.rows(), I.rows());         // default is column major
        poisson.reserve(Eigen::VectorXi::Constant(I.rows(), 7));

        for (int k = 0; k < I.rows(); k++) {
            for (int i = 0; i < 6; i++)
                poisson.insert(k, I(k, i)) = A(k, i);
            poisson.insert(k, k) = -1.0;
        }
        poisson.makeCompressed();

        Eigen::BiCGSTAB < Eigen::SparseMatrix < std::complex < double >> > solver(poisson);
        x = solver.solveWithGuess(b, x);
    }


};


class Label {
    int nw, nh;
    std::vector <std::string> labels;
public:
    Label() : nw(0), nh(0) {
        labels.clear();
    }

    inline int getNw() const {
        return nw;
    }

    inline int getNh() const {
        return nh;
    }

    inline bool isEmpty() const {
        return labels.size() == 0;
    }

    inline bool hasSameSize(const Label &label) const {
        return nw == label.nw && nh == label.nh && label.labels.size() == labels.size();
    }

    inline void clear() {
        nw = 0;
        nh = 0;
        labels.clear();
    }

    inline void getId(Eigen::VectorXi &id, const std::map<std::string, int> &idmap) const {

        id.resize(labels.size());
        for (int k = 0; k < labels.size(); k++) {
            auto it = idmap.find(labels[k]);
            if (it != idmap.end())
                id(k) = it->second;
            else id(k) = 0;
        }
    }

    inline void modifyById(const Eigen::VectorXi &id, const std::map<int, std::string> &idmap) {
        if (id.rows() != nw * nh) return;

        labels.resize(id.rows());
        for (int k = 0; k < labels.size(); k++) {
            auto it = idmap.find(id(k));
            if (it != idmap.end())
                labels[k] = it->second;
            else labels[k] = ".";

        }
    }

    inline void fromFile(const std::string &fn, const std::string &key) {

        clear();
        std::ifstream in(fn.c_str());
        fromStream(in, key);
        in.close();

        if (nw * nh != labels.size()) clear();

    }

    inline void fromString(const std::string &text, const std::string &key) {
        clear();
        std::istringstream in(text);
        fromStream(in, key);

        if (nw * nh != labels.size()) clear();

    }

    inline void fromStream(std::istream &is, const std::string &key) {
        clear();
        std::string findkey;
        while (findkey != key) if (!(is >> findkey)) break;

        if (findkey != key) {
            clear();
        }

        is >> nw >> nh;

        labels.resize(nw * nh);

        for (int j = nh - 1; j >= 0; j--)
            for (int i = 0; i < nw; i++) {
                int k = j * nw + i;
                if (!(is >> labels[k])) labels[k] = ".";
            }

        if (nw * nh != labels.size()) clear();

    }

    inline void toStream(std::ostream &os, const std::string &key) const {

        assert(nw * nh == labels.size());

        os << key << ' ' << nw << ' ' << nh << std::endl;

        for (int j = nh - 1; j >= 0; j--) {
            if (j % 2 != 0) os << ' ';
            for (int i = 0; i < nw; i++) {
                int k = j * nw + i;
                os << labels[k] << ' ';
            }
            os << std::endl;
        }


    }

    inline void toString(std::string &text, const std::string &key) const {

        std::ostringstream os;
        toStream(os, key);
        text = os.str();

    }

    inline void toFile(const std::string &fn, const std::string &key) const {

        std::ofstream os(fn);
        toStream(os, key);
        os.close();

    }
};

struct Boundary {
    double R;
    double T;
    std::vector<double> unode;
    std::vector<double> dtnode;
    std::vector<bool> where;
    std::string currentVar;

    Boundary() : R(0.0), currentVar("_") {
        unode.clear();
        dtnode.clear();
        where.clear();
    }

    inline double eval(const std::map<std::string, double> &variables) const {

        if (T == 0.0 || unode.size() == 0) return 0.0;
        if (unode.size() != dtnode.size()) throw std::runtime_error("Boundary node problem!");

        double t = 0.0;
        auto it = variables.find("t");
        if (it != variables.end()) t = it->second;
        t = std::fmod(t, T);

        double current = 0.0;
        it = variables.find(currentVar);
        if (it != variables.end()) current = it->second;

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

        return udrive - R * current;
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

template<typename T>
inline bool lineGetValueList(std::vector <T> &val, const std::vector <std::string> &lines,
                             const std::string commandkey, const std::vector <T> &defval) {
    val = defval;
    for (std::string line : lines) {
        std::istringstream ls(line);
        std::string command;

        ls >> command;
        if (command == commandkey) {
            int no;
            ls >> no;

            val.clear();
            for (int k = 0; k < no; k++) {
                T item;
                ls >> item;
                val.push_back(item);
            }
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
    idmap.clear();
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
ls >>
command;
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
ls >>
command;
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

    std::vector<std::tuple<std::vector < bool>, int, int, double, double, double>>      oxidation,
    reduction;
    std::vector <std::tuple<std::vector < std::tuple < int, char, int>>,int,int,double>>
    pattern_transition;

    std::vector <std::tuple<std::string, std::vector < bool>, std::vector<bool>>>
    current_surfaces;
    std::vector <std::tuple<std::string, std::vector < bool>>>
    voltage_sites;

    Parameters() {
        reset();
    }

    inline void reset() {

        a = 1.0;
        h = 1.0;
        dt = 1.0;
        kT = 1.0;
        omega = 0.001;
        lab2type["."] = 0;
        lab2site["."] = 0;
        type2lab[0] = ".";
        site2lab[0] = ".";
        tsource.resize(1);
        tsource.setZero();
        tcharge.resize(1);
        tcharge.setZero();
        ttau.resize(1);
        ttau.setOnes();
        tbarrier.resize(1, 1);
        tbarrier.setZero();
        tgamma1st.resize(1, 1);
        tgamma1st.setZero();
        tgamma2nd.resize(1, 1);
        tgamma2nd.setZero();
        ter.resize(1, 1);
        ter.setOnes();
        tsigma.resize(1, 1);
        tsigma.setZero();
        tec = std::complex<double>(1.0, 0.0) * ter - 4.0 * 3.1415 / omega * std::complex<double>(0.0, 1.0) * tsigma;
        oxidation.clear();
        reduction.clear();

        pattern_transition.clear();
        current_surfaces.clear();
        voltage_sites.clear();

        boundary.clear();
    }

    inline void logStream(std::ostream &os) const {
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

    inline void toStream(std::ostream &os) const {

        os << "dt " << dt << std::endl;
        os << "a " << a << std::endl;
        os << "h " << h << std::endl;
        os << "kT " << kT << std::endl;
        os << "omega " << omega << std::endl;


        os << "charge 0" << std::endl;
        os << "tau 0" << std::endl;
        os << "source 0" << std::endl;

        os << "er 0 " << ter << std::endl;
        os << "sigma 0 " << tsigma << std::endl;

        os << "gamma1st 0 " << tgamma1st << std::endl;
        os << "gamma2nd 0 " << tgamma2nd << std::endl;
        os << "barrier 0 " << tbarrier << std::endl;

        os << "type 1 . " << std::endl;
        os << "site 1 . " << std::endl;


        os << "oxidation 0" << std::endl;
        os << "reduction 0" << std::endl;
        os << "pattern 0" << std::endl;
        os << "current 0" << std::endl;
        os << "voltage 0" << std::endl;
        os << "boundary 0" << std::endl;

    }

    inline void toString(std::string &text) const {

        std::ostringstream os;
        toStream(os);
        text = os.str();

    }

    inline void toFile(const std::string &fn) const {

        std::ofstream os(fn);
        toStream(os);
        os.close();

    }


    inline void fromStream(std::istream &is) {

        std::vector <std::string> lines;
        std::string line;
        while (std::getline(is, line)) lines.push_back(line);


        if (!lineGetIdMap(lab2type, lines, "type")) {
            reset();
            return;
        }
        if (!lineGetIdMap(lab2site, lines, "site")) {
            reset();
            return;
        }

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

    inline void getIEfield(Eigen::Matrix<double, Eigen::Dynamic, 2> &IFi,
                           Eigen::Matrix<double, Eigen::Dynamic, 2> &EFi,
                           const Eigen::Matrix<double, Eigen::Dynamic, 6> &I,
                           const Eigen::Matrix<double, Eigen::Dynamic, 6> &E,
                           const Eigen::Matrix<double, 6, 2> &an, const Eigen::VectorXi &B) const {
        assert(I.rows() == E.rows());

        IFi.resize(I.rows(), 2);
        EFi.resize(I.rows(), 2);

        for (int k = 0; k < I.rows(); k++) {

            IFi(k, 0) = 0;
            IFi(k, 1) = 0;
            EFi(k, 0) = 0;
            EFi(k, 1) = 0;

            if (B(k) == 1) continue;

            for (int i = 0; i < 6; i++) {
                IFi(k, 0) += I(k, i) * an(i, 0);
                IFi(k, 1) += I(k, i) * an(i, 1);
                EFi(k, 0) += E(k, i) * an(i, 0);
                EFi(k, 1) += E(k, i) * an(i, 1);
            }

            IFi(k, 0) /= 6.0;
            IFi(k, 1) /= 6.0;
            EFi(k, 0) /= 6.0;
            EFi(k, 1) /= 6.0;
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
            double mured = efred - u(ari) - qs(ari) / dosred;
            double muoxMmured = muox - mured;
            for (int iox = 0; iox < 6; iox++)
                for (int ire = 0; ire < 6; ire++) {
                    int aoinb = nb_array(aoi, iox);
                    int arinb = nb_array(ari, ire);
                    muox = efox - u(aoinb) - qs(aoinb) / dosox;
                    mured = efred - u(arinb) - qs(arinb) / dosred;

                    muoxMmured = std::min(muoxMmured, muox - mured);
                }


            //double dEox = elchempot(aoi,taoito,typenb,u) - elchempot(aoi,taoifr,typenb,u) + muox;
            //double dEred = elchempot(ari,tarito,typenb,u) - elchempot(ari,tarifr,typenb,u) - mured;

            //double ee = dEox + dEred;

            double ee = elchempot(aoi, taoito, typenb, u) - elchempot(aoi, taoifr, typenb, u) +
                        elchempot(ari, tarito, typenb, u) - elchempot(ari, tarifr, typenb, u) +
                        muoxMmured;

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

    Eigen::VectorXi types;
    Eigen::VectorXi sites;
    int nw, nh;

    Eigen::VectorXi B;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> u0;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> u;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 6> e;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> b;

    Eigen::Matrix<double, Eigen::Dynamic, 6> I;
    Eigen::Matrix<double, Eigen::Dynamic, 6> E;
    Eigen::Matrix<double, Eigen::Dynamic, 1> q;

    Eigen::Matrix<double, Eigen::Dynamic, 2> IFi;
    Eigen::Matrix<double, Eigen::Dynamic, 2> EFi;

    MultiGridPoissonSolver mgp;
    SparseMatrixSolver smp;

    bool flagPoisson;
    bool flagDiffusion;
    bool flagRedox;
    bool flagPattern;

    bool flagDoOnlyBICG;

    std::vector <std::string> color_state;
    std::vector <std::string> color_sites;

    Control() {
        initVariables();
        generateGrid(nb, co, an, 0, 0);
        nw = 0;
        nh = 0;
    }

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

        assert(state.hasSameSize(site));

        nw = state.getNw();
        nh = state.getNh();

        if (nw == 0 || nh == 0) return;
        assert(nw > 0 && nh > 0);

        generateGrid(nb, co, an, state.getNw(), state.getNh());

        state.getId(types, parameters.lab2type);
        site.getId(sites, parameters.lab2site);

        parameters.getboundary(B, u0, sites, variables);
        u = u0;

        parameters.geteb(e, b, nb, types, sites, variables);

        //sparseMatrixSolver(u,nb,e,b);
        if (flagDoOnlyBICG) smp.solveFirst(u, nb, e, b); else mgp.solveFirst(u, e, -b, nw, nh);

        parameters.getIEq(I, E, q, u, nb, types);

        parameters.getVarCurrent(variables, I, nb, sites);

        parameters.getVarVoltage(variables, u.real(), sites);
        //std::cout << "--" << std::endl;
    }

    inline void step(const Parameters &parameters) {
        if (flagDiffusion) parameters.diffusion(types, u.real(), nb, generator);
        if (flagRedox) parameters.redox(types, u.real(), q, nb, generator);
        if (flagPattern) parameters.pattern(types, nb, generator);
        if (flagPoisson) {
            parameters.geteb(e, b, nb, types, sites, variables);
            //sparseMatrixSolver(u,nb,e,b);
            //std::cout << u(41) <<std::endl;

            if (flagDoOnlyBICG) smp.solve(u, nb, e, b); else mgp.solve(u, e, -b);
            std::cout << u(41) << std::endl;
            std::cout << "----------" << std::endl;
            parameters.getIEq(I, E, q, u, nb, types);
            parameters.getVarCurrent(variables, I, nb, sites);
            parameters.getVarVoltage(variables, u.real(), sites);
        }

        variables["t"] += parameters.dt;
    }

    inline void logStream(std::ostream &os, const Parameters &parameters) {

        os << variables["t"] << '\t';
        for (auto v : variables) if (v.first != "_" && v.first != "t") os << v.second << '\t';
        os << "var" << std::endl;
        //os << nw << ' ' << nh << std::endl;
        os << co.col(0).transpose() * parameters.a << std::endl;
        os << co.col(1).transpose() * parameters.a << std::endl;
        os << types.transpose() << std::endl;
        os << parameters.type2lab.size() << std::endl;
        for (int t = 0; t < parameters.type2lab.size(); t++) {
            os << color_state[t] << ' ';
            for (int k = 0; k < types.rows(); k++) if (types(k) == t) os << k << ' ';
            os << std::endl;
        }
        os << 6 << std::endl;
        os << u.real().transpose() << std::endl;
        os << q.transpose() << std::endl;

        parameters.getIEfield(IFi, EFi, I, E, an, B);

        os << IFi.col(0).transpose() << std::endl;
        os << IFi.col(1).transpose() << std::endl;
        os << EFi.col(0).transpose() << std::endl;
        os << EFi.col(1).transpose() << std::endl;

        os << std::endl;
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

        lineGetValueList(color_state, lines, "colorState",
                         {"lightgray", "black", "red", "gray", "darkblue", "darkred", "lightgray", "black", "red",
                          "gray", "darkblue", "darkred"});
        lineGetValueList(color_sites, lines, "colorSites",
                         {"lightgray", "black", "red", "gray", "darkblue", "darkred", "lightgray", "black", "red",
                          "gray", "darkblue", "darkred"});

        lineGetValue(mgp.MG_FULL_NCYCLES_FINE, lines, "poisson_maxiter", 50);
        lineGetValue(mgp.MG_FINAL_TOLERANCE, lines, "poisson_tolerance", 1.0e-5);
        lineGetValue(mgp.doBICG, lines, "poisson_doBICG", true);

        lineGetValue(flagDoOnlyBICG, lines, "poisson_doOnlyBICG", false);

        double time;
        lineGetValue(time, lines, "time", 0.0);
        variables["t"] = time;

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

#endif // GRWOTH_BASE_H

/*
    void prolong(Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> &fu,
                 const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> &u, int nw,int nh) {
        fu.resize(nw*nh*4);
        for(int j = 0; j < nh; j++)
            for(int i = 0; i < nw; i++) {
                int k = index(i,j,nw,nh);
                fu(index(2*i,2*j,2*nw,2*nh)) = u(index(i,j,nw,nh));
                fu(index(2*i+1,2*j,2*nw,2*nh)) = 0.5 * ( u(index(i,j,nw,nh)) + u(index(i+1,j,nw,nh)) );
                fu(index(2*i,2*j+1,2*nw,2*nh)) = 0.5 * (u(index(i,j,nw,nh)) + u(index(i,j+1,nw,nh)) );
                fu(index(2*i+1,2*j+1,2*nw,2*nh)) = 0.5 * (u(index(i+1,j,nw,nh)) + u(index(i,j+1,nw,nh)) );
            }
    }

    void generate_ladder(std::vector<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1>> & b_ladder,
                         std::vector<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1>> & e_ladder) {

    }

    void FMG(Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> &u,
             const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,6> &e,
             const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> &b,
             const Eigen::Matrix<int,Eigen::Dynamic,6> & nb6,int nw,int nh) {

        if(nw%2 != 0 || nh%2 != 0) {
            exactSolver(u,e,b,nb6);
            //for(int iter = 0; iter<100; iter++) smooth(u,e,b,nb6,nw,nh);
            return;
        }
        smooth(u,e,b,nb6,nw,nh);

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> r;
        residual(r,u,e,b,nb6);

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> bc;
        restrict(bc,r,nw,nh);

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,6> ec;
        restrict(ec,e,nb6,nw,nh);

        Eigen::Matrix<int,Eigen::Dynamic,6> nb6c;
        gen_nb6(nb6c,nw/2,nh/2);

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> vc(nb6c.rows()); vc.setZero();
        FMG(vc,ec,bc,nb6c,nw/2,nh/2);

        Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> v;
        prolong(v,vc,nw/2,nh/2);
        u += v;
        smooth(u,e,b,nb6,nw,nh);

    }
template<typename T>
    void exactSolver(Eigen::Matrix<T,Eigen::Dynamic,1> & u,
                     const Eigen::Matrix<T,Eigen::Dynamic,6> & e,
                     const Eigen::Matrix<T,Eigen::Dynamic,1> & b,
                     const Eigen::Matrix<int,Eigen::Dynamic,6> & nb6) {
        if(b.rows() == 0 || e.rows() == 0 || nb6.rows() == 0) return;

        if(u.rows() != b.rows()) {
            u.resize(b.rows());
            u.setZero();
        }

        Eigen::SparseMatrix<T> poisson(nb6.rows(),nb6.rows());         // default is column major
        poisson.reserve(Eigen::VectorXi::Constant(nb6.rows(),7));


        for(int k = 0; k<nb6.rows(); k++ ) {
            for(int i = 0; i < 6; i++)
                poisson.insert(k,nb6(k,i)) = e(k,i);
            poisson.insert(k,k) = -1.0;
        }
        poisson.makeCompressed();

        Eigen::BiCGSTAB<Eigen::SparseMatrix<T>> solver(poisson);

        u = solver.solveWithGuess(-b,u);

    }

    void residual(Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> &r,
                  const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> &u,
                  const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,6> &e,
                  const Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> &b,
                  const Eigen::Matrix<int,Eigen::Dynamic,6> & nb6) {
        r.resizeLike(u);
        for(int k = 0; k < nb6.rows(); k++) {
            std::complex<double> sum = b(k) - u(k);
            for(int l = 0; l < 6; l++) sum += u(nb6(k,l))*e(k,l);
            r(k) = sum;
        }
    }

*/
