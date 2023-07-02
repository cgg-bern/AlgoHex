/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

// implementation based on David Palmer's ARFF implementation, MultiSdp.cpp

#include <AlgoHex/Config/Version.hh>

#if ALGOHEX_WITH_MOSEK

#include "SHProjectorSDP.hh"
#include "SDPMatrices.hh"
#include "OctaToQuat.hh"

#include <fusion.h>


using mosek::fusion::Model;
using mosek::fusion::Domain;
using mosek::fusion::Expr;
using mosek::fusion::Set;
using mosek::fusion::ObjectiveSense;
using monty::shape;
using monty::ndarray;
using monty::new_array_ptr;

namespace AlgoHex {


template<typename EigenT>
std::shared_ptr<ndarray<double, 2>> fromEigen(const EigenT &m)
{
    auto rows = m.rows();
    auto cols = m.cols();
    auto A = new_array_ptr<double, 2>(shape(rows, cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            (*A)(i, j) = m(i,j);
        }
    }
    return A;
}

struct SHProjectorSDP::impl {
    impl()
        : model(new Model())
    {}
    ~impl() {model->dispose();}

    Model::t model;
};

AlgoHex::SHProjectorSDP::SHProjectorSDP()
    : impl_(std::make_unique<impl>())
{
    static auto A = fromEigen(sdpA);
    static auto b = fromEigen(sdpB);

    // we prefer multithreading in thr outer loop;
    model()->setSolverParam("numThreads", 1);
    model()->setSolverParam("intpntMultiThread", "off");
    model()->setSolverParam("intpntCoTolRelGap", "1.0e-12");

    const auto d = 9;
    auto Q = model()->variable("Q", Domain::inPSDCone(d + 1));
    model()->constraint(Expr::mul(A, Q->reshape(Set::make((d + 1) * (d + 1), 1))), Domain::equalsTo(b));

}

AlgoHex::SHProjectorSDP::~SHProjectorSDP()
{
}

SHProjectionResult SHProjectorSDP::project(const SHCoeffs &sh_coeffs)
{
    SHProjectionResult res;
    res.shc = project_only_shc(sh_coeffs);
    res.q = from_projected_shc(res.shc);
    return res;
}

AlgoHex::SHCoeffs AlgoHex::SHProjectorSDP::project_only_shc(const AlgoHex::SHCoeffs &shc)
{
    const auto d = 9;

    auto Q = model()->getVariable("Q");

    auto q0 = new_array_ptr<double, 2>(shape(d, 1));
    for (int i = 0; i < d; ++i) {
        (*q0)[i] = shc[i];
    }
    auto qi = Q->slice(Set::make(1, 0), Set::make(d + 1, 1));

    model()->objective(ObjectiveSense::Minimize,
            Expr::add(Expr::sum(Q->diag()),
                      Expr::mul(-2.0, Expr::dot(q0, qi))));

    model()->solve();
    auto qValue = *(qi->level());
    SHCoeffs result;
    for (int i = 0; i < d; ++i) {
        result[i] = qValue[i];
    }
    // TODO: check if this was successful?
    return result;
}

mosek::fusion::Model *SHProjectorSDP::model()
{
    return impl_->model.get();
}

} // namespace AlgoHex

#endif // ALGOHEX_WITH_MOSEK
