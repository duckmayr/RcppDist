// RcppDist.h
//
// Copyright (C) 2018 JB Duck-Mayr
//
// This file is part of RcppDist.
//
// RcppDist is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppDist is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppDist.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RCPPDIST_RCPPDIST_H
#define RCPPDIST_RCPPDIST_H

#ifdef RCPPDIST_DONT_USE_ARMA
    #include <Rcpp.h>
#else
    #include <RcppArmadillo.h>
    #include <mvnorm.h>
    #include <mvt.h>
    #include <wishart.h>
#endif

#include <4beta.h>
#include <lst.h>
#include <truncnorm.h>
#include <trunct.h>
#include <trunclst.h>
#include <triangular.h>

#endif
