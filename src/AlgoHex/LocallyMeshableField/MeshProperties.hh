/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include <AlgoHex/TypeDef.hh>

namespace AlgoHex
{
template<class MeshT>
class MeshPropertiesT
{
public:
  MeshPropertiesT(MeshT &_mesh) :
          mesh_(_mesh),
          valence_(mesh_.template request_edge_property<int>("edge_valance")),
          trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
          cell_quaternions_(mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions")),
          sgl_vt_(mesh_.template request_vertex_property<int>("singular_vertex")),
          target_length_(mesh_.template request_vertex_property<double>("target_length")),
          feature_node_(mesh_.template request_vertex_property<int>("AlgoHex::FeatureVertices")),
          feature_edge_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureEdgeVertices")),
          feature_edge_(mesh_.template request_edge_property<int>("AlgoHex::FeatureEdges")),
          feature_face_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices")),
          feature_face_edge_(mesh_.template request_edge_property<bool>("AlgoHex::FeatureFaceEdges")),
          feature_fprop_(mesh_.template request_face_property<int>("AlgoHex::FeatureFaces"))
  {
  }


  void initialize_singular_vertex_property()
  {
    mesh_.set_persistent(sgl_vt_, true);

    for (const auto &vh: mesh_.vertices())
      sgl_vt_[vh] = 0;

    for (const auto &ehi: mesh_.edges())
    {
      if (valence_[ehi] != 0)
      {
        if (!mesh_.is_boundary(ehi))
        {
          sgl_vt_[mesh_.edge(ehi).from_vertex()] = 1;
          sgl_vt_[mesh_.edge(ehi).to_vertex()] = 1;
        }
      }
    }
    for (const auto &ehi: mesh_.edges())
    {
      if (valence_[ehi] != 0)
      {
        if (mesh_.is_boundary(ehi))
        {
          sgl_vt_[mesh_.edge(ehi).from_vertex()] = 2;
          sgl_vt_[mesh_.edge(ehi).to_vertex()] = 2;
        }
      }
    }
  }

//        void import_feature_properties() {
//            auto vcolor = mesh_.template request_vertex_property<int>("vertex_colors");
//            auto ecolor = mesh_.template request_edge_property<int>("edge_colors");
//            auto fcolor = mesh_.template request_face_property<int>("face_colors");
//
//            mesh_.set_persistent(feature_node_, true);
//            mesh_.set_persistent(feature_edge_, true);
//            mesh_.set_persistent(feature_fprop_, true);
//
//            for(const auto vhi : mesh_.vertices()) {
//                if(vcolor[vhi] != 0)
//                    feature_node_[vhi] = true;
//                else
//                    feature_node_[vhi] = false;
//            }
//
//            for(const auto ehi : mesh_.edges()) {
//                if(ecolor[ehi] != 0)
//                    feature_edge_[ehi] = true;
//                else
//                    feature_edge_[ehi] = false;
//            }
//
//            for(const auto fhi : mesh_.faces()) {
//                if(fcolor[fhi] != 0)
//                    feature_fprop_[fhi] = true;
//                else
//                    feature_fprop_[fhi] = false;
//            }
//        }
  void import_feature_properties()
  {
    auto vcolor = mesh_.template request_vertex_property<int>("vertex_colors");
    auto ecolor = mesh_.template request_edge_property<int>("edge_colors");
    auto fcolor = mesh_.template request_face_property<int>("face_colors");

    mesh_.set_persistent(feature_node_, true);
    mesh_.set_persistent(feature_edge_, true);
    mesh_.set_persistent(feature_fprop_, true);

    for (const auto vhi: mesh_.vertices())
      feature_node_[vhi] = vcolor[vhi];

    for (const auto ehi: mesh_.edges())
      feature_edge_[ehi] = ecolor[ehi];

    for (const auto fhi: mesh_.faces())
      feature_fprop_[fhi] = fcolor[fhi];
  }

  void initialize_feature_vertex_property()
  {
    mesh_.set_persistent(feature_edge_vertex_, true);
    mesh_.set_persistent(feature_face_edge_, true);
    mesh_.set_persistent(feature_face_vertex_, true);

    //reset
    for (const auto vhi: mesh_.vertices())
    {
      feature_edge_vertex_[vhi] = false;
      feature_face_vertex_[vhi] = false;
    }
    for (const auto ehi: mesh_.edges())
    {
      feature_face_edge_[ehi] = false;
    }

    //assign
    for (const auto &ehi: mesh_.edges())
    {
      if (feature_edge_[ehi] > 0)
      {
        feature_edge_vertex_[mesh_.edge(ehi).from_vertex()] = true;
        feature_edge_vertex_[mesh_.edge(ehi).to_vertex()] = true;
      }
    }

    //update feature node property if not prescribed
    auto it_me = std::max_element(feature_node_.begin(), feature_node_.end());
    if (it_me != feature_node_.end() && *it_me == 0)
    {
      for (const auto &vhi: mesh_.vertices())
      {
        if (feature_edge_vertex_[vhi])
        {
          if (n_incident_feature_edges(vhi) > 2)
            feature_node_[vhi] = 1;
        }
      }
    }


    for (const auto &fhi: mesh_.faces())
    {
      if (feature_fprop_[fhi] > 0)
      {
        for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
        {
          feature_face_edge_[*fe_it] = true;

          feature_face_vertex_[mesh_.edge(*fe_it).from_vertex()] = true;
          feature_face_vertex_[mesh_.edge(*fe_it).to_vertex()] = true;
        }
      }
    }
  }

  void initialize_target_edge_length()
  {
    mesh_.set_persistent(target_length_, true);

    for (const auto vh: mesh_.vertices())
    {
      double len(0.);
      int n = 0;
      for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
      {
        len += mesh_.length(*ve_it);
        n++;
      }
      len /= (double) n;

      target_length_[vh] = len;
    }
  }

  void smooth_target_edge_length(const int _iter = 5, const double _damping = 0.8)
  {
    for (int i = 0; i < _iter; ++i)
    {
      for (const auto &vh: mesh_.vertices())
      {
        if (sgl_vt_[vh] == 0)
        {
          double len = 0.;
          int n = 0;
          for (auto vv_it = mesh_.vv_iter(vh); vv_it.valid(); ++vv_it)
          {
            len += target_length_[*vv_it];
            n++;
          }
          len /= (double) n;

          target_length_[vh] = len;
        }
        else
        {
          double len = 0.;
          int n = 0;
          for (auto vv_it = mesh_.vv_iter(vh); vv_it.valid(); ++vv_it)
          {
            len += target_length_[*vv_it];
            n++;
          }

          len /= (double) n;

          target_length_[vh] = _damping * target_length_[vh] + (1 - _damping) * len;
        }
      }
    }
  }

  void localize_properties()
  {
    mesh_.set_persistent(sgl_vt_, false);

    mesh_.set_persistent(target_length_, false);
  }

  int n_incident_feature_edges(const VH _vh) const
  {
    int n_fe = 0;
    for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
      if (feature_edge_[*ve_it] > 0)
        n_fe++;

    return n_fe;
  }

  void update_singular_vertex_property(const VH _vh)
  {
    int n_b_sge = 0, n_i_sge = 0;
    for (auto voe_it = mesh_.ve_iter(_vh); voe_it.valid(); ++voe_it)
      if (valence_[*voe_it] != 0)
      {
        if (mesh_.is_boundary(*voe_it))
          n_b_sge++;
        else
          n_i_sge++;
      }

    if (n_i_sge > 0)
      sgl_vt_[_vh] = 1;
    if (n_b_sge > 0)
      sgl_vt_[_vh] = 2;
    if (n_i_sge == 0 && n_b_sge == 0)
      sgl_vt_[_vh] = 0;
  }


  //heuristically get the node index by counting edge valence
  int node_index(const VH _vh) const
  {
    if (!mesh_.is_boundary(_vh))
    {
      int n_all = 0, n_val_ng1_i = 0, n_val1_i = 0, n_val_ng2_i = 0, n_val2_i = 0, n_val3_i = 0, n_val4_i = 0, n_complex = 0;
      for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
      {
        auto eh = mesh_.edge_handle(*voh_it);
        if (valence_[eh] == -1)
          n_val_ng1_i++;
        else if (valence_[eh] == 1)
          n_val1_i++;
        else if (valence_[eh] == -2)
          n_val_ng2_i++;
        else if (valence_[eh] == 2)
          n_val2_i++;
        else if (valence_[eh] == 3)
          n_val3_i++;
        else if (valence_[eh] == 4)
          n_val4_i++;
        else if (valence_[eh] != 0)
          n_complex++;

        if (valence_[eh] != 0)
          n_all++;
      }

      if (n_complex > 0)
        return 20;

      if (n_val_ng2_i != 0 || n_val2_i != 0 || n_val3_i != 0 || n_val4_i != 0)
      {
        if (n_val_ng1_i == 0 && n_val1_i == 0 && n_val_ng2_i == 0 && n_val2_i == 2)
          return 0;
        else if (n_val_ng1_i == 0 && n_val1_i == 0 && n_val_ng2_i == 2 && n_val2_i == 0)
          return 0;
        else if (n_all == 2 && n_val3_i == 2)
          return 0;
        else if (n_all == 2 && n_val4_i == 2)
          return 0;
        else if (n_all == 2 && n_val_ng2_i == 1 && n_val2_i == 1)//zipper node of valence +2 and -2
          return 12;
        else
          return 20;
      }

      if ((n_val_ng1_i == 0 && n_val1_i == 0) || (n_val_ng1_i == 2 && n_val1_i == 0) ||
          (n_val_ng1_i == 0 && n_val1_i == 2))
        return 0;
      else if (n_val_ng1_i == 4 && n_val1_i == 0)
        return 1;
      else if (n_val_ng1_i == 2 && n_val1_i == 2)
        return 2;
      else if (n_val_ng1_i == 1 && n_val1_i == 3)
        return 3;
      else if (n_val_ng1_i == 0 && n_val1_i == 4)
        return 4;
      else if (n_val_ng1_i == 2 && n_val1_i == 6)
        return 5;
      else if (n_val_ng1_i == 0 && n_val1_i == 6)
        return 6;
      else if (n_val_ng1_i == 0 && n_val1_i == 8)
        return 7;
      else if (n_val_ng1_i == 0 && n_val1_i == 12)
        return 8;
      else if (n_val_ng1_i == 1 && n_val1_i == 1)
      {
        return 10;
      }
      else
        return 20;
    }
    else
    {
      int n_all = 0, n_val_ng1_b = 0, n_val_ng2_b = 0, n_val1_b = 0, n_val2_b = 0,
              n_val2_i = 0, n_val3_i = 0, n_val4_i = 0, n_val_ng2_i = 0, n_val_ng1_i = 0, n_val1_i = 0, n_complex = 0;
      for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
      {
        EH eh = mesh_.edge_handle(*voh_it);
        if (mesh_.is_boundary(eh))
        {
          if (valence_[eh] == -1)
            n_val_ng1_b++;
          else if (valence_[eh] == 1)
            n_val1_b++;
          else if (valence_[eh] == -2)
            n_val_ng2_b++;
          else if (valence_[eh] == 2)
            n_val2_b++;
        }
        else
        {
          if (valence_[eh] == -1)
            n_val_ng1_i++;
          else if (valence_[eh] == 1)
            n_val1_i++;
          else if (valence_[eh] == -2)
            n_val_ng2_i++;
          else if (valence_[eh] == 2)
            n_val2_i++;
          else if (valence_[eh] == 3)
            n_val3_i++;
          else if (valence_[eh] == 4)
            n_val4_i++;
          else if (valence_[eh] != 0)
            n_complex++;
        }

        if (valence_[eh] != 0)
          n_all++;
      }

      if (n_complex > 0)
        return 20;

      if (n_val_ng2_i != 0 || n_val2_i != 0 || n_val3_i != 0 || n_val4_i != 0)
      {
        if (n_val_ng1_b == 0 && n_val1_b == 0 && n_val_ng1_i == 0 && n_val1_i == 0 &&
            n_val_ng2_b == 0 && n_val2_b == 0 && n_val_ng2_i == 1 && n_val2_i == 0)
          return -10;
        else if (n_val_ng1_b == 0 && n_val1_b == 0 && n_val_ng1_i == 0 && n_val1_i == 0 &&
                 n_val_ng2_b == 0 && n_val2_b == 0 && n_val_ng2_i == 0 && n_val2_i == 1)
          return -11;
        else if (n_all == 1 && n_val3_i == 1) // only enumerate simple case
          return -18;
        else if (n_all == 1 && n_val4_i == 1) // only enumerate simple case
          return -19;
        else
          return 20;
      }

      //all cases that num_sge <=4
      if (n_val_ng2_b == 0 && n_val2_b == 0)
      {
        if (((n_val_ng1_b == 0 && n_val1_b == 0 && n_val_ng1_i == 0 && n_val1_i == 0) ||
             (n_val_ng1_b == 2 && n_val1_b == 0 && n_val_ng1_i == 0 && n_val1_i == 0) ||
             (n_val_ng1_b == 0 && n_val1_b == 2 && n_val_ng1_i == 0 && n_val1_i == 0)))
          return 0;
        else if ((n_val_ng1_b == 3 && n_val1_b == 0) && n_val_ng1_i == 0 && n_val1_i == 0)
          return -1;
        else if ((n_val_ng1_b == 2 && n_val1_b == 1) && n_val_ng1_i == 0 && n_val1_i == 0)
          return -2;
        else if ((n_val_ng1_b == 2 && n_val1_b == 2) && n_val_ng1_i == 0 && n_val1_i == 0)
          return -3; // -4 has the same number, but it's a mirror case
        else if ((n_val_ng1_b == 0 && n_val1_b == 3) && n_val_ng1_i == 0 && n_val1_i == 1)
          return -5;
        else if ((n_val_ng1_b == 1 && n_val1_b == 2) && n_val_ng1_i == 0 && n_val1_i == 0)
          return -6;
        else if ((n_val_ng1_b == 0 && n_val1_b == 3) && n_val_ng1_i == 0 && n_val1_i == 0)
          return -7;
        else if ((n_val_ng1_b == 1 && n_val1_b == 2) && n_val_ng1_i == 1 && n_val1_i == 0)
          return -8;
          //this could be an invalid case if the interior arc goes in the wrong direction
        else if ((n_val_ng1_b == 0 && n_val1_b == 2) && n_val_ng1_i == 0 && n_val1_i == 1)
          return -9;
        else if (n_val_ng1_b == 0 && n_val1_b == 0 && n_val_ng1_i == 1 && n_val1_i == 0)
          return -10;
        else if (n_val_ng1_b == 0 && n_val1_b == 0 && n_val_ng1_i == 0 && n_val1_i == 1)
          return -11;
          //this could be an invalid case if the interior arc goes in the wrong direction
        else if (n_val_ng1_b == 0 && n_val1_b == 2 && n_val_ng1_i == 1 && n_val1_i == 0)
          return -12;
        else if (n_val_ng1_b == 0 && n_val1_b == 2 && n_val_ng1_i == 1 && n_val1_i == 1)
          return -13;
        else if ((n_val_ng1_b == 1 && n_val1_b == 2) && n_val_ng1_i == 0 && n_val1_i == 1)
          return -14;
        else if (((n_val_ng1_b == 1 && n_val1_b == 1) && n_val_ng1_i == 0 && n_val1_i == 0) ||
                 ((n_val_ng1_b == 1 && n_val1_i == 1) && n_val_ng1_i == 0 && n_val1_b == 0) ||
                 ((n_val_ng1_i == 1 && n_val1_b == 1) && n_val_ng1_b == 0 && n_val1_i == 0))
          return 10;
        else
          return 20;
      }
      else
      {//valence +2
        if (n_val_ng1_i == 0 && n_val1_i == 0 && n_val_ng1_b == 0 && n_val1_b == 0)
        {
          if ((n_val_ng2_b == 2 && n_val2_b == 0) ||
              (n_val_ng2_b == 0 && n_val2_b == 2))
          {
            return 0;
          }
        }

        //two 0-sector touch
        if (n_val_ng2_b == 0 && n_val2_b == 1 && n_val_ng1_b == 2 && n_val1_b == 2 && n_val_ng1_i == 0 && n_val1_i == 0)
          return -15;

        if (n_val_ng2_b == 0 && n_val2_b == 1 && n_val_ng1_b == 2 && n_val1_b == 0 && n_val_ng1_i == 0 && n_val1_i == 0)
          return -16;

        if (n_val_ng2_b == 0 && n_val2_b == 1 && n_val_ng1_b == 0 && n_val1_b == 2 && n_val_ng1_i == 0 && n_val1_i == 0)
          return -17;

        return 20;
      }
    }
  }

  bool is_singular_node(const VH _vh) const
  {
    int n = 0;
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
      if (valence_[mesh_.edge_handle(*voh_it)] != 0)
        n++;

    if (n > 2)
      return true;
    return false;
  }


protected:
  MeshT &mesh_;

  //input
  EP<int> valence_;
  HFP<int> trans_prop_;
  CP <Quaternion> cell_quaternions_;

  VP<int> sgl_vt_;
  VP<double> target_length_;

  VP<int> feature_node_;

  VP<bool> feature_edge_vertex_;
  EP<int> feature_edge_;

  VP<bool> feature_face_vertex_;
  EP<bool> feature_face_edge_;
  FP<int> feature_fprop_;
};
}



