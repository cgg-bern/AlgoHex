/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <AlgoHex/TypeDef.hh>
#include "EnergyFunctions.hh"

namespace AlgoHex
{
class SingularVertexOptProblemQuaternion : public COMISO::NProblemInterface
{
public:
  using Edge = std::pair<int, int>;
  using Face = std::array<int, 3>;
  using VTuple = std::array<int, 3>;
  using Cell = std::array<int, 4>;
  using DualEdge = std::pair<int, int>;

  using Mat = Mat3d;

  SingularVertexOptProblemQuaternion(const std::vector<double> &_points, std::vector<double> &_partial_points) :
          points_(_points), partial_points_(_partial_points)
  {
    x_imrme_.resize(12);

    g_imrme2_.resize(6);
    h_imrme2_.resize(6, 6);

    g_imrme3_.resize(9);
    h_imrme3_.resize(9, 9);

    g_imrme4_.resize(12);
    h_imrme4_.resize(12, 12);

    x_simrme_.resize(9);

    x_nme_.resize(12);

    x_pe_.resize(PerpendicularEnergyCPRD::n_unknowns());
    g_pe_.resize(PerpendicularEnergyCPRD::n_unknowns());
    h_pe_.resize(PerpendicularEnergyCPRD::n_unknowns(), PerpendicularEnergyCPRD::n_unknowns());

    x_cpe_.resize(ComplexEdgeLength::n_unknowns());
    g_cpe_.resize(ComplexEdgeLength::n_unknowns());
    h_cpe_.resize(ComplexEdgeLength::n_unknowns(), ComplexEdgeLength::n_unknowns());

    x_cve_.resize(CurvatureSmooth::n_unknowns());
    g_cve_.resize(CurvatureSmooth::n_unknowns());
    h_cve_.resize(CurvatureSmooth::n_unknowns(), CurvatureSmooth::n_unknowns());

    x_rpe_.resize(RepulsionEnergy::n_unknowns());
  }

  ~SingularVertexOptProblemQuaternion() {}

  //
  virtual int n_unknowns() { return partial_points_.size(); }

  virtual void initial_x(double *_x)
  {
    for (auto i = 0u; i < partial_points_.size(); ++i)
      _x[i] = partial_points_[i];
  }

  double initial_f()
  {
    std::vector<double> x(n_unknowns());
    initial_x(x.data());
    return eval_f(x.data());
  }

  virtual double eval_f(const double *_x)
  {
    double energy(0);

    double ea = energy;
    //perpendicular energy
    for (auto j = 0u; j < ppd_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_pe_[i] = _x[3 * ppd_edges_[j].first + i];
        x_pe_[i + 3] = _x[3 * ppd_edges_[j].second + i];
      }

      double ee = ppd_weights_[j] * PerpendicularEnergyCPRD::eval_f(ppd_axes_[j], x_pe_.data());

//                if(!std::isfinite(ee))
//                    std::cout << "Error: alignmennt e infinite! "<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<" w: "<<ppd_weights_[j]<< " f: "
//                              <<ee <<" c pt: "<<x_pe_.transpose() << " vhs: "<<ppd_edges_[j].first<<" "<<ppd_edges_[j].second;

      energy += ee;
    }

//            std::cout << " alignmennt energy: " << energy - ea;
//            ea = energy;

    //complex edge length energy
    for (auto j = 0u; j < comp_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cpe_[i] = _x[3 * comp_edges_[j].first + i];
        x_cpe_[i + 3] = _x[3 * comp_edges_[j].second + i];
      }
      double ee = cp_weights_[j] * ComplexEdgeLength::eval_f(x_cpe_.data());

//                if(!std::isfinite(ee))
//                    std::cout << "Error: complex e infinite! "<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<" w: "<<cp_weights_[j]<< " f: "
//                              <<ee <<" c pt: "<<x_cpe_.transpose() << " vhs: "<<comp_edges_[j].first<<" "<<comp_edges_[j].second;

      energy += ee;
    }

//            std::cout << " complex edge length energy: " << energy - ea;
//            ea = energy;

    //curvature smooth energy
    for (auto j = 0u; j < vertex_tuples_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cve_[i] = _x[3 * vertex_tuples_[j][0] + i];
        x_cve_[i + 3] = _x[3 * vertex_tuples_[j][1] + i];
        x_cve_[i + 6] = _x[3 * vertex_tuples_[j][2] + i];
      }

      double ee = cv_weights_[j] * CurvatureSmooth::eval_f(x_cve_.data());

//                if(!std::isfinite(ee))
//                    std::cout << "Error: curvature smooth e infinite! "<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<" w: "<<cv_weights_[j]<< " f: "
//                          <<ee <<" c pt: "<<x_cve_.transpose() << " vhs: "<<vertex_tuples_[j][0]<<" "<<vertex_tuples_[j][1]<<" "<<vertex_tuples_[j][2];

      energy += ee;
    }

//            std::cout << " curvature smoothing energy: " << energy - ea;
//            ea = energy;

    //repulsion smooth energy
    for (auto j = 0u; j < rp_vertices_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_rpe_[i] = _x[3 * rp_vertices_[j] + i];
      }

      energy += rp_weights_[j] * RepulsionEnergy::eval_f(x_rpe_.data(), target_points_[j]);
    }

//            std::cout << " repulsion energy: " << energy - ea;
//            ea = energy;

    //inverse mean ratio metric
    int w = 0;
    for (const auto &c: cells1_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_imrme_[j] = _x[3 * c[0] + j];

      //take from the original points
      for (int i = 1; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      //
      double ee = imrm_weights1_[w] * IMRMEnergy::eval_f(x_imrme_.data());

//                if(!std::isfinite(ee))
//                    std::cout << "Error: imrm1_ e infinite! "<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<" w: "<<imrm_weights1_[w]<< " f: "
//                              <<ee <<" c pt: "<<x_imrme_.transpose() << " vhs: "<<c[0]<<" "<<c[1]<<" "<<c[2]<<" "<<c[3];


      w++;
      energy += ee;
    }

//            std::cout<<" cell1 vol energy: "<<energy-ea;
//            ea = energy;

    w = 0;
    for (const auto &c: cells2_)
    {
      //take from the unknown
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 2; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      //
      double e = imrm_weights2_[w] * IMRMEnergy2::eval_f(x_imrme_.data());

      w++;

      energy += e;
    }

//            std::cout<<" e cell2 vol energy: "<<energy-ea;
//            ea = energy;

    w = 0;
    for (const auto &c: cells3_)
    {
      //take from the unknown
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 3; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      //
      double e = imrm_weights3_[w] * IMRMEnergy3::eval_f(x_imrme_.data());

      w++;

      energy += e;
    }

//            std::cout<<" cell3 vol energy: "<<energy-ea;
//            ea = energy;


    w = 0;
    for (const auto &c: cells4_)
    {
      //take from the unknown
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //
      double e = imrm_weights4_[w] * IMRMEnergy4::eval_f(x_imrme_.data());

      w++;

      energy += e;
    }
//            std::cout << " cell 4 vol energy: " << energy - ea;
//            ea = energy;


    //surface inverse mean ratio metric
    w = 0;
    for (const auto &f: faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_simrme_[j] = _x[3 * f[0] + j];

      //take from the original points
      Vec3d pp, px, p0(x_simrme_[0], x_simrme_[1], x_simrme_[2]);
      for (int i = 1; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
          px[j] = points_[3 * f[i] + j];

        project_on_tangent_plane(p0, px, normals_[w], pp);

        for (int j = 0; j < 3; ++j)
        {
          x_simrme_[3 * i + j] = pp[j];
        }
      }
      //
      double e = simrm_weights_[w] * SIMRMEnergy::eval_f(x_simrme_.data());

      w++;
      energy += e;
    }

//            std::cout << " surface imrm energy: " << energy - ea;
//            ea = energy;

    //create a tetrahedron on top of the triangle
    //if the normal of the face flips, the tet volume will be negative
    w = 0;
    for (const auto &f: base_faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_nme_[j] = _x[3 * f[0] + j];

      //take from the original points
      for (int i = 1; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_nme_[3 * i + j] = points_[3 * f[i] + j];
        }

      //the fourth point
      for (int j = 0; j < 3; ++j)
        x_nme_[9 + j] = top_pts_[w][j];

      double ee = c_weights_[w] * IMRMEnergy::eval_f(x_nme_.data());

//                if(!std::isfinite(ee))
//                    std::cout << "Error: normal flip e infinite! "<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<" w: "<<normal_weights_[w]<< " f: "
//                              <<ee <<" face pt: "<<x_nme_.transpose();

      energy += ee;
      w++;
    }
//            std::cout << " anti normal  energy: " << energy - ea << std::endl;


    return energy;
  }


  double initial_f(const double *_x)
  {
    double energy(0);

    double ea = energy;
    //perpendicular energy
    for (auto j = 0u; j < ppd_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_pe_[i] = _x[3 * ppd_edges_[j].first + i];
        x_pe_[i + 3] = _x[3 * ppd_edges_[j].second + i];
      }

      double ee = ppd_weights_[j] * PerpendicularEnergyCPRD::eval_f(ppd_axes_[j], x_pe_.data());

      if (!std::isfinite(ee))
        std::cout << "Error: alignmennt e infinite! "
                  << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << " w: " << ppd_weights_[j]
                  << " f: "
                  << ee << " c pt: " << x_pe_.transpose() << " vhs: " << ppd_edges_[j].first << " "
                  << ppd_edges_[j].second;

      energy += ee;
    }

    std::cout << " alignmennt energy: " << energy - ea;
    ea = energy;

    //complex edge length energy
    for (auto j = 0u; j < comp_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cpe_[i] = _x[3 * comp_edges_[j].first + i];
        x_cpe_[i + 3] = _x[3 * comp_edges_[j].second + i];
      }
      double ee = cp_weights_[j] * ComplexEdgeLength::eval_f(x_cpe_.data());

      if (!std::isfinite(ee))
        std::cout << "Error: edge shrinking e infinite! "
                  << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << " w: " << cp_weights_[j]
                  << " f: "
                  << ee << " c pt: " << x_cpe_.transpose() << " vhs: " << comp_edges_[j].first << " "
                  << comp_edges_[j].second;

      energy += ee;
    }

    std::cout << " edge shrinking energy: " << energy - ea;
    ea = energy;

    //curvature smooth energy
    for (auto j = 0u; j < vertex_tuples_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cve_[i] = _x[3 * vertex_tuples_[j][0] + i];
        x_cve_[i + 3] = _x[3 * vertex_tuples_[j][1] + i];
        x_cve_[i + 6] = _x[3 * vertex_tuples_[j][2] + i];
      }

      double ee = cv_weights_[j] * CurvatureSmooth::eval_f(x_cve_.data());

      if (!std::isfinite(ee))
        std::cout << "Error: curvature smooth e infinite! "
                  << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << " w: " << cv_weights_[j]
                  << " f: "
                  << ee << " c pt: " << x_cve_.transpose() << " vhs: " << vertex_tuples_[j][0] << " "
                  << vertex_tuples_[j][1] << " " << vertex_tuples_[j][2];

      energy += ee;
    }

    std::cout << " curvature smoothing energy: " << energy - ea;
    ea = energy;

    //repulsion smooth energy
    for (auto j = 0u; j < rp_vertices_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_rpe_[i] = _x[3 * rp_vertices_[j] + i];
      }

      energy += rp_weights_[j] * RepulsionEnergy::eval_f(x_rpe_.data(), target_points_[j]);
    }

    std::cout << " repulsion energy: " << energy - ea;
    ea = energy;

    //inverse mean ratio metric
    int w = 0;
    for (const auto &c: cells1_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_imrme_[j] = _x[3 * c[0] + j];

      //take from the original points
      for (int i = 1; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      //
      double ee = imrm_weights1_[w] * IMRMEnergy::eval_f(x_imrme_.data());

      if (!std::isfinite(ee))
        std::cout << "Error: imrm1_ e infinite! " << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
                  << " w: " << imrm_weights1_[w] << " f: "
                  << ee << " c pt: " << x_imrme_.transpose() << " vhs: " << c[0] << " " << c[1] << " " << c[2] << " "
                  << c[3] << std::endl;


      w++;
      energy += ee;
    }

    std::cout << " cell1 vol energy: " << energy - ea;
    ea = energy;

    w = 0;
    for (const auto &c: cells2_)
    {
      //take from the unknown
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 2; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      //
      double ee = imrm_weights2_[w] * IMRMEnergy2::eval_f(x_imrme_.data());

      if (!std::isfinite(ee))
        std::cout << "Error: imrm2_ e infinite! " << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
                  << " w: " << imrm_weights2_[w] << " f: "
                  << ee << " c pt: " << x_imrme_.transpose() << " vhs: " << c[0] << " " << c[1] << " " << c[2] << " "
                  << c[3] << std::endl;

      w++;

      energy += ee;
    }

    std::cout << " e cell2 vol energy: " << energy - ea;
    ea = energy;

    w = 0;
    for (const auto &c: cells3_)
    {
      //take from the unknown
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 3; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      //
      double ee = imrm_weights3_[w] * IMRMEnergy3::eval_f(x_imrme_.data());
      if (!std::isfinite(ee))
      {
        std::cout << "Error: imrm3_ e infinite! "
                  << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << " w: "
                  << imrm_weights3_[w] << " f: "
                  << ee << " c pt: " << x_imrme_.transpose() << " vhs: " << c[0] << " " << c[1] << " "
                  << c[2] << " " << c[3];

        {
          auto v01 = Vec3d(x_imrme_[3] - x_imrme_[0], x_imrme_[4] - x_imrme_[1], x_imrme_[5] - x_imrme_[2]);
          auto v02 = Vec3d(x_imrme_[6] - x_imrme_[0], x_imrme_[7] - x_imrme_[1], x_imrme_[8] - x_imrme_[2]);
          auto v03 = Vec3d(x_imrme_[9] - x_imrme_[0], x_imrme_[10] - x_imrme_[1], x_imrme_[11] - x_imrme_[2]);

          double vprod = x_imrme_[0] * x_imrme_[5] * x_imrme_[7] - x_imrme_[0] * x_imrme_[4] * x_imrme_[8] +
                         x_imrme_[1] * x_imrme_[3] * x_imrme_[8] - x_imrme_[1] * x_imrme_[5] * x_imrme_[6] -
                         x_imrme_[2] * x_imrme_[3] * x_imrme_[7] + x_imrme_[2] * x_imrme_[4] * x_imrme_[6] +
                         x_imrme_[0] * x_imrme_[4] * x_imrme_[11] - x_imrme_[0] * x_imrme_[5] * x_imrme_[10] -
                         x_imrme_[1] * x_imrme_[3] * x_imrme_[11] + x_imrme_[1] * x_imrme_[5] * x_imrme_[9] +
                         x_imrme_[2] * x_imrme_[3] * x_imrme_[10] - x_imrme_[2] * x_imrme_[4] * x_imrme_[9] -
                         x_imrme_[0] * x_imrme_[7] * x_imrme_[11] + x_imrme_[0] * x_imrme_[8] * x_imrme_[10] +
                         x_imrme_[1] * x_imrme_[6] * x_imrme_[11] - x_imrme_[1] * x_imrme_[8] * x_imrme_[9] -
                         x_imrme_[2] * x_imrme_[6] * x_imrme_[10] + x_imrme_[2] * x_imrme_[7] * x_imrme_[9] +
                         x_imrme_[3] * x_imrme_[7] * x_imrme_[11] - x_imrme_[3] * x_imrme_[8] * x_imrme_[10] -
                         x_imrme_[4] * x_imrme_[6] * x_imrme_[11] + x_imrme_[4] * x_imrme_[8] * x_imrme_[9] +
                         x_imrme_[5] * x_imrme_[6] * x_imrme_[10] - x_imrme_[5] * x_imrme_[7] * x_imrme_[9];
          std::cout << " eignenvol: " << v01.cross(v02).dot(v03) << " vol: " << vprod;


          double cdet = -x_imrme_[0] * x_imrme_[5] * x_imrme_[10] + x_imrme_[0] * x_imrme_[8] * x_imrme_[10] +
                        x_imrme_[0] * x_imrme_[4] * x_imrme_[11] - x_imrme_[0] * x_imrme_[7] * x_imrme_[11] -
                        x_imrme_[0] * x_imrme_[4] * x_imrme_[8] + x_imrme_[0] * x_imrme_[5] * x_imrme_[7] -
                        x_imrme_[1] * x_imrme_[3] * x_imrme_[11] + x_imrme_[1] * x_imrme_[6] * x_imrme_[11] +
                        x_imrme_[1] * x_imrme_[3] * x_imrme_[8] - x_imrme_[1] * x_imrme_[5] * x_imrme_[6] +
                        x_imrme_[1] * x_imrme_[5] * x_imrme_[9] - x_imrme_[1] * x_imrme_[8] * x_imrme_[9] +
                        x_imrme_[2] * x_imrme_[3] * x_imrme_[10] - x_imrme_[2] * x_imrme_[6] * x_imrme_[10] -
                        x_imrme_[3] * x_imrme_[8] * x_imrme_[10] + x_imrme_[5] * x_imrme_[6] * x_imrme_[10] +
                        x_imrme_[3] * x_imrme_[7] * x_imrme_[11] - x_imrme_[4] * x_imrme_[6] * x_imrme_[11] -
                        x_imrme_[2] * x_imrme_[3] * x_imrme_[7] + x_imrme_[2] * x_imrme_[4] * x_imrme_[6] -
                        x_imrme_[2] * x_imrme_[4] * x_imrme_[9] + x_imrme_[2] * x_imrme_[7] * x_imrme_[9] +
                        x_imrme_[4] * x_imrme_[8] * x_imrme_[9] - x_imrme_[5] * x_imrme_[7] * x_imrme_[9];
          std::cout << " det: " << cdet;
        }
        std::cout << std::endl;
      }

      w++;

      energy += ee;
    }

    std::cout << " cell3 vol energy: " << energy - ea;
    ea = energy;


    w = 0;
    for (const auto &c: cells4_)
    {
      //take from the unknown
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //
      double ee = imrm_weights4_[w] * IMRMEnergy4::eval_f(x_imrme_.data());


      if (!std::isfinite(ee))
        std::cout << "Error: imrm4_ e infinite! " << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
                  << " w: " << imrm_weights4_[w] << " f: "
                  << ee << " c pt: " << x_imrme_.transpose() << " vhs: " << c[0] << " " << c[1] << " " << c[2] << " "
                  << c[3] << std::endl;

      w++;

      energy += ee;
    }
    std::cout << " cell 4 vol energy: " << energy - ea;
    ea = energy;


    //surface inverse mean ratio metric
    w = 0;
    for (const auto &f: faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_simrme_[j] = _x[3 * f[0] + j];

      //take from the original points
      Vec3d pp, px, p0(x_simrme_[0], x_simrme_[1], x_simrme_[2]);
      for (int i = 1; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
          px[j] = points_[3 * f[i] + j];

        project_on_tangent_plane(p0, px, normals_[w], pp);

        for (int j = 0; j < 3; ++j)
        {
          x_simrme_[3 * i + j] = pp[j];
        }
      }
      //
      double ee = simrm_weights_[w] * SIMRMEnergy::eval_f(x_simrme_.data());

      if (!std::isfinite(ee))
        std::cout << "Error: simrm e infinite! " << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
                  << " w: " << simrm_weights_[w] << " f: "
                  << ee << " c pt: " << x_simrme_.transpose() << " vhs: " << f[0] << " " << f[1] << " " << f[2]
                  << std::endl;

      w++;
      energy += ee;
    }

    std::cout << " surface imrm energy: " << energy - ea;
    ea = energy;


    //create a tetrahedron on top of the triangle
    //if the normal of the face flips, the tet volume will be negative
    w = 0;
    for (const auto &f: base_faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_nme_[j] = _x[3 * f[0] + j];

      //take from the original points
      for (int i = 1; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_nme_[3 * i + j] = points_[3 * f[i] + j];
        }

      //the fourth point
      for (int j = 0; j < 3; ++j)
        x_nme_[9 + j] = top_pts_[w][j];

      double ee = c_weights_[w] * IMRMEnergy::eval_f(x_nme_.data());

//                if(!std::isfinite(ee))
//                    std::cout << "Error: normal flip e infinite! "<<std::setprecision(std::numeric_limits<long double>::digits10 + 1)<<" w: "<<normal_weights_[w]<< " f: "
//                              <<ee <<" face pt: "<<x_nme_.transpose();

      energy += ee;
      w++;
    }
    std::cout << " anti normal flip energy: " << energy - ea << std::endl;


    return energy;
  }


  //project px to the plane defined by p0 and n
  template<typename Vec3>
  static void project_on_tangent_plane(const Vec3 &_p0, const Vec3 &_px, const Vec3 &_n,
                                       Vec3 &_pp)
  {
    _pp = _px - (_px - _p0).dot(_n) * _n;
  }


  virtual void eval_gradient(const double *_x, double *_g)
  {
    // set to zero
    memset(_g, 0, this->n_unknowns() * sizeof(double));

    //perpendicular energy
    for (auto j = 0u; j < ppd_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_pe_[i] = _x[3 * ppd_edges_[j].first + i];
        x_pe_[i + 3] = _x[3 * ppd_edges_[j].second + i];
      }

      PerpendicularEnergyCPRD::eval_gradient(ppd_axes_[j], x_pe_.data(), g_pe_.data());

      //copy to global
      for (unsigned int ii = 0; ii < 3; ++ii)
      {
        _g[3 * ppd_edges_[j].first + ii] += ppd_weights_[j] * g_pe_[ii];
        _g[3 * ppd_edges_[j].second + ii] += ppd_weights_[j] * g_pe_[3 + ii];
        if (!std::isfinite(ppd_weights_[j]) || !std::isfinite(g_pe_[ii]) || !std::isfinite(g_pe_[ii + 3]))
        {
          std::cout << "Error: ppd gradient infinite!" << std::endl;
        }
      }
    }


    //complex edge length
    for (auto j = 0u; j < comp_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cpe_[i] = _x[3 * comp_edges_[j].first + i];
        x_cpe_[i + 3] = _x[3 * comp_edges_[j].second + i];
      }

      ComplexEdgeLength::eval_gradient(x_cpe_.data(), g_cpe_.data());

      //copy to global
      for (unsigned int ii = 0; ii < 3; ++ii)
      {
        _g[3 * comp_edges_[j].first + ii] += cp_weights_[j] * g_cpe_[ii];
        _g[3 * comp_edges_[j].second + ii] += cp_weights_[j] * g_cpe_[3 + ii];
      }
    }


    //curvature smooth energy
    for (auto j = 0u; j < vertex_tuples_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cve_[i] = _x[3 * vertex_tuples_[j][0] + i];
        x_cve_[i + 3] = _x[3 * vertex_tuples_[j][1] + i];
        x_cve_[i + 6] = _x[3 * vertex_tuples_[j][2] + i];
      }

      CurvatureSmooth::eval_gradient(x_cve_.data(), g_cve_.data());

      //copy to global
      for (unsigned int ii = 0; ii < 3; ++ii)
      {
        _g[3 * vertex_tuples_[j][0] + ii] += cv_weights_[j] * g_cve_[ii];
        _g[3 * vertex_tuples_[j][1] + ii] += cv_weights_[j] * g_cve_[3 + ii];
        _g[3 * vertex_tuples_[j][2] + ii] += cv_weights_[j] * g_cve_[6 + ii];
      }
    }

    //repulsion smooth energy
    for (auto j = 0u; j < rp_vertices_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_rpe_[i] = _x[3 * rp_vertices_[j] + i];
      }

      RepulsionEnergy::eval_gradient(x_rpe_.data(), target_points_[j], g_rpe_.data());

      for (unsigned int ii = 0; ii < 3; ++ii)
        _g[3 * rp_vertices_[j] + ii] += rp_weights_[j] * g_rpe_[ii];
    }


    //inverse mean ratio metric
    int w = 0;
    for (const auto &c: cells1_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_imrme_[j] = _x[3 * c[0] + j];

      //take from the original points
      for (int i = 1; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      IMRMEnergy::eval_gradient(x_imrme_.data(), g_imrme1_);

      //copy to global
      for (int j = 0; j < 3; ++j)
      {
        double gtmp = imrm_weights1_[w] * g_imrme1_[j];
        if (!std::isfinite(gtmp))
          gtmp = 0;

        _g[3 * c[0] + j] += gtmp;
      }

      w++;
    }

    w = 0;
    for (const auto &c: cells2_)
    {
      //take from the unknown
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 2; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      IMRMEnergy2::eval_gradient(x_imrme_.data(), g_imrme2_.data());

      //copy to global
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
        {
          double gtmp = imrm_weights2_[w] * g_imrme2_[3 * i + j];
          if (!std::isfinite(gtmp))
            gtmp = 0;
          _g[3 * c[i] + j] += gtmp;
        }

      w++;
    }


    w = 0;
    for (const auto &c: cells3_)
    {
      //take from the unknown
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 3; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      IMRMEnergy3::eval_gradient(x_imrme_.data(), g_imrme3_.data());

      //copy to global
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          double gtmp = imrm_weights3_[w] * g_imrme3_[3 * i + j];
          if (!std::isfinite(gtmp))
            gtmp = 0;

          _g[3 * c[i] + j] += gtmp;
        }

      w++;
    }

    w = 0;
    for (const auto &c: cells4_)
    {
      //take from the unknown
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      IMRMEnergy4::eval_gradient(x_imrme_.data(), g_imrme4_.data());

      //copy to global
//                bool err = false;
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          double gtmp = imrm_weights4_[w] * g_imrme4_[3 * i + j];
          if (!std::isfinite(gtmp))
          {
            gtmp = 0;
//                            err = true;
          }
          _g[3 * c[i] + j] += gtmp;
        }

      w++;
    }


    w = 0;
    for (const auto &f: faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_simrme_[j] = _x[3 * f[0] + j];

      //take from the original points
      Vec3d pp, px, p0(x_simrme_[0], x_simrme_[1], x_simrme_[2]);
      for (int i = 1; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
          px[j] = points_[3 * f[i] + j];

        project_on_tangent_plane(p0, px, normals_[w], pp);

        for (int j = 0; j < 3; ++j)
        {
          x_simrme_[3 * i + j] = pp[j];
        }
      }
      //
      SIMRMEnergy::eval_gradient(x_simrme_.data(), g_simrme_);


      //copy to global
      for (int j = 0; j < 3; ++j)
      {
        double gtmp = simrm_weights_[w] * g_simrme_[j];
        if (!std::isfinite(gtmp))
          gtmp = 0;

        _g[3 * f[0] + j] += gtmp;
      }

      w++;
    }

    //create a tetrahedron on top of the triangle
    //if the normal of the face flips, the tet volume will be negative
    w = 0;
    for (const auto &f: base_faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_nme_[j] = _x[3 * f[0] + j];

      //take from the original points
      for (int i = 1; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_nme_[3 * i + j] = points_[3 * f[i] + j];
        }

      //the fourth point
      for (int j = 0; j < 3; ++j)
        x_nme_[9 + j] = top_pts_[w][j];

      IMRMEnergy::eval_gradient(x_nme_.data(), g_nme_);

      //copy to global
      for (int j = 0; j < 3; ++j)
      {
        double gtmp = c_weights_[w] * g_nme_[j];
        if (!std::isfinite(gtmp))
          gtmp = 0;

        _g[3 * f[0] + j] += gtmp;
      }

      w++;
    }

  }

  virtual void eval_hessian(const double *_x, SMatrixNP &_H)
  {
    _H.resize(n_unknowns(), n_unknowns());
    _H.setZero();

    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(
            36 * ppd_edges_.size() + 36 * comp_edges_.size() + 9 * cells1_.size() + 36 * cells2_.size() +
            +81 * cells3_.size() + +324 * cells4_.size() + +9 * faces_.size() + 3 * rp_vertices_.size() +
            9 * base_faces_.size());


    //perpendicular energy
    for (auto j = 0u; j < ppd_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_pe_[i] = _x[3 * ppd_edges_[j].first + i];
        x_pe_[i + 3] = _x[3 * ppd_edges_[j].second + i];
      }

      PerpendicularEnergyCPRD::eval_hessian(ppd_axes_[j], x_pe_.data(), h_pe_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
        {
          triplets.emplace_back(3 * ppd_edges_[j].first + m, 3 * ppd_edges_[j].first + n,
                                ppd_weights_[j] * h_pe_(m, n));
          triplets.emplace_back(3 * ppd_edges_[j].second + m, 3 * ppd_edges_[j].second + n,
                                ppd_weights_[j] * h_pe_(m + 3, n + 3));
          triplets.emplace_back(3 * ppd_edges_[j].first + m, 3 * ppd_edges_[j].second + n,
                                ppd_weights_[j] * h_pe_(m, n + 3));
          triplets.emplace_back(3 * ppd_edges_[j].second + m, 3 * ppd_edges_[j].first + n,
                                ppd_weights_[j] * h_pe_(m + 3, n));
        }
    }

    //complex edge energy
    for (auto j = 0u; j < comp_edges_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cpe_[i] = _x[3 * comp_edges_[j].first + i];
        x_cpe_[i + 3] = _x[3 * comp_edges_[j].second + i];
      }

      ComplexEdgeLength::eval_hessian(x_cpe_.data(), h_cpe_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
        {
          triplets.emplace_back(3 * comp_edges_[j].first + m, 3 * comp_edges_[j].first + n,
                                cp_weights_[j] * h_cpe_(m, n));
          triplets.emplace_back(3 * comp_edges_[j].second + m, 3 * comp_edges_[j].second + n,
                                cp_weights_[j] * h_cpe_(m + 3, n + 3));
          triplets.emplace_back(3 * comp_edges_[j].first + m, 3 * comp_edges_[j].second + n,
                                cp_weights_[j] * h_cpe_(m, n + 3));
          triplets.emplace_back(3 * comp_edges_[j].second + m, 3 * comp_edges_[j].first + n,
                                cp_weights_[j] * h_cpe_(m + 3, n));
        }
    }


    //curvature smooth energy
    for (auto j = 0u; j < vertex_tuples_.size(); ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        x_cve_[i] = _x[3 * vertex_tuples_[j][0] + i];
        x_cve_[i + 3] = _x[3 * vertex_tuples_[j][1] + i];
        x_cve_[i + 6] = _x[3 * vertex_tuples_[j][2] + i];
      }

      CurvatureSmooth::eval_hessian(x_cve_.data(), h_cve_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
        {
          triplets.emplace_back(3 * vertex_tuples_[j][m], 3 * vertex_tuples_[j][n],
                                cv_weights_[j] * h_cve_(3 * m, 3 * n));
          triplets.emplace_back(3 * vertex_tuples_[j][m] + 1, 3 * vertex_tuples_[j][n] + 1,
                                cv_weights_[j] * h_cve_(3 * m + 1, 3 * n + 1));
          triplets.emplace_back(3 * vertex_tuples_[j][m] + 2, 3 * vertex_tuples_[j][n] + 2,
                                cv_weights_[j] * h_cve_(3 * m + 2, 3 * n + 2));
        }
    }


    //repulsion smooth energy
    for (auto j = 0u; j < rp_vertices_.size(); ++j)
    {
      RepulsionEnergy::eval_hessian(h_rpe_);

      for (unsigned int ii = 0; ii < 3; ++ii)
        triplets.emplace_back(3 * rp_vertices_[j] + ii, 3 * rp_vertices_[j] + ii,
                              rp_weights_[j] * h_rpe_(ii, ii));
    }


    //inverse mean ratio metric
    int w = 0;
    for (const auto &c: cells1_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_imrme_[j] = _x[3 * c[0] + j];

      //take from the original points
      for (int i = 1; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      IMRMEnergy::eval_hessian(x_imrme_.data(), h_imrme1_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          triplets.emplace_back(3 * c[0] + m, 3 * c[0] + n, imrm_weights1_[w] * h_imrme1_(m, n));

      w++;
    }

    w = 0;
    for (const auto &c: cells2_)
    {
      //take from the unknown
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 2; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      IMRMEnergy2::eval_hessian(x_imrme_.data(), h_imrme2_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
              triplets.emplace_back(3 * c[i] + m, 3 * c[j] + n,
                                    imrm_weights2_[w] * h_imrme2_(m + 3 * i, n + 3 * j));

      w++;
    }

    w = 0;
    for (const auto &c: cells3_)
    {
      //take from the unknown
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      //take from the original points
      for (int i = 3; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_imrme_[3 * i + j] = points_[3 * c[i] + j];
        }

      IMRMEnergy3::eval_hessian(x_imrme_.data(), h_imrme3_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
              triplets.emplace_back(3 * c[i] + m, 3 * c[j] + n,
                                    imrm_weights3_[w] * h_imrme3_(m + 3 * i, n + 3 * j));

      w++;
    }


    w = 0;
    for (const auto &c: cells4_)
    {
      //take from the unknown
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
          x_imrme_[3 * i + j] = _x[3 * c[i] + j];

      IMRMEnergy4::eval_hessian(x_imrme_.data(), h_imrme4_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
              triplets.emplace_back(3 * c[i] + m, 3 * c[j] + n,
                                    imrm_weights4_[w] * h_imrme4_(m + 3 * i, n + 3 * j));

      w++;
    }



    //surface imrm
    w = 0;
    for (const auto &f: faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_simrme_[j] = _x[3 * f[0] + j];

      //take from the original points
      Vec3d pp, px, p0(x_simrme_[0], x_simrme_[1], x_simrme_[2]);
      for (int i = 1; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
          px[j] = points_[3 * f[i] + j];

        project_on_tangent_plane(p0, px, normals_[w], pp);

        for (int j = 0; j < 3; ++j)
        {
          x_simrme_[3 * i + j] = pp[j];
        }
      }

      SIMRMEnergy::eval_hessian(x_simrme_.data(), h_simrme_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          triplets.emplace_back(3 * f[0] + m, 3 * f[0] + n, simrm_weights_[w] * h_simrme_(m, n));

      w++;
    }


    //create a tetrahedron on top of the triangle
    //if the normal of the face flips, the tet volume will be negative
    w = 0;
    for (const auto &f: base_faces_)
    {
      //take from the unknown
      for (int j = 0; j < 3; ++j)
        x_nme_[j] = _x[3 * f[0] + j];

      //take from the original points
      for (int i = 1; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          x_nme_[3 * i + j] = points_[3 * f[i] + j];
        }

      //the fourth point
      for (int j = 0; j < 3; ++j)
        x_nme_[9 + j] = top_pts_[w][j];

      IMRMEnergy::eval_hessian(x_nme_.data(), h_nme_);

      //copy to global
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          triplets.emplace_back(3 * f[0] + m, 3 * f[0] + n, c_weights_[w] * h_nme_(m, n));

      w++;
    }

    _H.setFromTriplets(triplets.begin(), triplets.end());
  }

  virtual void store_result(const double *_x)
  {
    for (auto i = 0u; i < partial_points_.size(); ++i)
      partial_points_[i] = _x[i];
  }

  //one unknown vertex
  void add_imrm_energy_element_one(const int _v_newid0, const int _vid1, const int _vid2, const int _vid3,
                                   const double _weight)
  {
    cells1_.push_back(Cell{_v_newid0, _vid1, _vid2, _vid3});
    imrm_weights1_.push_back(_weight);
  }

  //two
  void add_imrm_energy_element_two(const int _v_newid0, const int _v_newid1, const int _vid2, const int _vid3,
                                   const double _weight)
  {
    cells2_.push_back(Cell{_v_newid0, _v_newid1, _vid2, _vid3});
    imrm_weights2_.push_back(_weight);
  }

  //three
  void
  add_imrm_energy_element_three(const int _v_newid0, const int _v_newid1, const int _v_newid2, const int _vid3,
                                const double _weight)
  {
    cells3_.push_back(Cell{_v_newid0, _v_newid1, _v_newid2, _vid3});
    imrm_weights3_.push_back(_weight);
  }

  //four
  void
  add_imrm_energy_element_four(const int _v_newid0, const int _v_newid1, const int _v_newid2, const int _v_newid3,
                               const double _weight)
  {
    cells4_.push_back(Cell{_v_newid0, _v_newid1, _v_newid2, _v_newid3});
    imrm_weights4_.push_back(_weight);
  }

  void
  add_boundary_simrm_element(const int _v_newid0, const int _vid1, const int _vid2, const Vec3d &_n,
                             const double _weight)
  {
    faces_.push_back(Face{_v_newid0, _vid1, _vid2});
    normals_.push_back(_n);
    simrm_weights_.push_back(_weight);
  }

  void add_anti_normal_flip_element(const int _v_newid0, const int _vid1, const int _vid2, const Vec3d &_top_pt,
                                    const double _weight)
  {
    base_faces_.push_back(Face{_v_newid0, _vid1, _vid2});
    top_pts_.push_back(_top_pt);
    c_weights_.push_back(_weight);
  }


  void add_perpendicular_element(const int _vid0, const int _vid1, const Vec3d &_n, const double _w)
  {
    ppd_edges_.emplace_back(_vid0, _vid1);
    ppd_axes_.push_back(_n);
    ppd_weights_.push_back(_w);
  }

  void add_shrink_edge_element(const int _vid0, const int _vid1, const double _w)
  {
    comp_edges_.emplace_back(_vid0, _vid1);
    cp_weights_.push_back(_w);
  }

  void add_curvature_smooth_element(const int _vid0, const int _vid1, const int _vid2, const double _w)
  {
    vertex_tuples_.emplace_back(VTuple{_vid0, _vid1, _vid2});
    cv_weights_.push_back(_w);
  }

  void add_repulsion_term(const int _vid, const Vec3d &_target_point, const double _w)
  {
    rp_vertices_.push_back(_vid);
    target_points_.push_back(_target_point);
    rp_weights_.push_back(_w);
  }


private:
  //complete mesh vertex position
  const std::vector<double> &points_;
  //partial mesh vertex position
  std::vector<double> &partial_points_;

  //inverse mean ratio metric
  std::vector<Cell> cells1_, cells2_, cells3_, cells4_;
  std::vector<double> imrm_weights1_, imrm_weights2_, imrm_weights3_, imrm_weights4_;

  //
  Eigen::VectorXd x_imrme_;

  Vec3d g_imrme1_;
  Mat3d h_imrme1_;

  Eigen::VectorXd g_imrme2_;
  Eigen::MatrixXd h_imrme2_;

  Eigen::VectorXd g_imrme3_;
  Eigen::MatrixXd h_imrme3_;

  Eigen::VectorXd g_imrme4_;
  Eigen::MatrixXd h_imrme4_;


  //surface inverse mean ratio
  std::vector<Face> faces_;
  std::vector<Vec3d> normals_;
  std::vector<double> simrm_weights_;
  //
  Eigen::VectorXd x_simrme_;
  Vec3d g_simrme_;
  Mat3d h_simrme_;


  //anti normal flip energy
  std::vector<Face> base_faces_;
  std::vector<Vec3d> top_pts_;
  std::vector<double> c_weights_;
  //
  Eigen::VectorXd x_nme_;
  Vec3d g_nme_;
  Mat3d h_nme_;


  //axis perpendicular to edge direction (dot or cross product)
  std::vector<Edge> ppd_edges_;
  std::vector<Vec3d> ppd_axes_;
  std::vector<double> ppd_weights_;

  Eigen::VectorXd x_pe_;
  Eigen::VectorXd g_pe_;
  Eigen::MatrixXd h_pe_;

  //edge length term
  std::vector<Edge> comp_edges_;
  std::vector<double> cp_weights_;

  Eigen::VectorXd x_cpe_;
  Eigen::VectorXd g_cpe_;
  Eigen::MatrixXd h_cpe_;


  //curvature smoothing
  std::vector<VTuple> vertex_tuples_;
  std::vector<double> cv_weights_;

  Eigen::VectorXd x_cve_;
  Eigen::VectorXd g_cve_;
  Eigen::MatrixXd h_cve_;


  //repulsion term
  std::vector<int> rp_vertices_;
  std::vector<Vec3d> target_points_;
  std::vector<double> rp_weights_;

  Eigen::VectorXd x_rpe_;
  Vec3d g_rpe_;
  Mat3d h_rpe_;
};

}






