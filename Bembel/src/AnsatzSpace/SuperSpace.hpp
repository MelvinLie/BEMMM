// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_INCLUDE_SUPERSPACE_H_
#define BEMBEL_INCLUDE_SUPERSPACE_H_
namespace Bembel {
/**
 *  \ingroup AnsatzSpace
 *  \brief The superspace manages local polynomial bases on each element of the
 * mesh and provides an itnerface to evaluate them.
 */
template <typename Derived>
struct SuperSpace {
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  SuperSpace(){};
  SuperSpace(Geometry& geom, int M, int P) { init_SuperSpace(geom, M, P); }
  SuperSpace(const SuperSpace& other) {
    mesh_ = other.mesh_;
    phi = other.phi;
    phiDx = other.phiDx;
    phiPhi = other.phiPhi;
    phiPhiDx = other.phiPhiDx;
    phiPhiDy = other.phiPhiDy;
    phiTimesPhi = other.phiTimesPhi;
    // vPhiScalVPhi = other.vPhiScalVPhi;
    divPhiTimesDivPhi = other.divPhiTimesDivPhi;
    curlPhiTimesCurlPhi = other.curlPhiTimesCurlPhi;
    vecTimesCurlPhi = other.vecTimesCurlPhi;
    vecTimesCurlPhi_single_basis = other.vecTimesCurlPhi_single_basis;
    vecTimesCurlPhi_BSpline = other.vecTimesCurlPhi_BSpline;
    polynomial_degree = other.polynomial_degree;
    polynomial_degree_plus_one_squared =
        other.polynomial_degree_plus_one_squared;
  }
  SuperSpace(SuperSpace&& other) {
    mesh_ = other.mesh_;
    phi = other.phi;
    phiDx = other.phiDx;
    phiPhi = other.phiPhi;
    phiPhiDx = other.phiPhiDx;
    phiPhiDy = other.phiPhiDy;
    phiTimesPhi = other.phiTimesPhi;
    // vPhiScalVPhi = other.vPhiScalVPhi;
    divPhiTimesDivPhi = other.divPhiTimesDivPhi;
    curlPhiTimesCurlPhi = other.curlPhiTimesCurlPhi;
    vecTimesCurlPhi = other.vecTimesCurlPhi;
    vecTimesCurlPhi_single_basis = other.vecTimesCurlPhi_single_basis;
    vecTimesCurlPhi_BSpline = other.vecTimesCurlPhi_BSpline;
    polynomial_degree = other.polynomial_degree;
    polynomial_degree_plus_one_squared =
        other.polynomial_degree_plus_one_squared;
  };
  SuperSpace& operator=(SuperSpace other) {
    mesh_ = other.mesh_;
    phi = other.phi;
    phiDx = other.phiDx;
    phiPhi = other.phiPhi;
    phiPhiDx = other.phiPhiDx;
    phiPhiDy = other.phiPhiDy;
    phiTimesPhi = other.phiTimesPhi;
    // vPhiScalVPhi = other.vPhiScalVPhi;
    divPhiTimesDivPhi = other.divPhiTimesDivPhi;
    curlPhiTimesCurlPhi = other.curlPhiTimesCurlPhi;
    vecTimesCurlPhi = other.vecTimesCurlPhi;
    vecTimesCurlPhi_single_basis = other.vecTimesCurlPhi_single_basis;
    vecTimesCurlPhi_BSpline = other.vecTimesCurlPhi_BSpline;
    polynomial_degree = other.polynomial_degree;
    polynomial_degree_plus_one_squared =
        other.polynomial_degree_plus_one_squared;
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  inline int get_polynomial_degree() const { return polynomial_degree; }
  inline int get_polynomial_degree_plus_one_squared() const {
    return polynomial_degree_plus_one_squared;
  }
  inline int get_refinement_level() const { return mesh_->get_max_level(); }
  inline int get_number_of_elements() const {
    return mesh_->get_number_of_elements();
  }
  inline int get_number_of_patches() const {
    return mesh_->get_geometry().size();
  }
  const PatchVector& get_geometry() const { return mesh_->get_geometry(); }
  inline const ClusterTree& get_mesh() const { return *mesh_; };
  //////////////////////////////////////////////////////////////////////////////
  //    init_SuperSpace
  //////////////////////////////////////////////////////////////////////////////
  void init_SuperSpace(const Geometry& geom, int M, int P) {
    polynomial_degree = P;
    polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);
    phi = (Basis::BasisHandler<
           typename LinearOperatorTraits<Derived>::Scalar>::funPtrPhi(P));
    phiDx = (Basis::BasisHandler<
             typename LinearOperatorTraits<Derived>::Scalar>::funPtrPhiDx(P));
    phiPhi = (Basis::BasisHandler<
              typename LinearOperatorTraits<Derived>::Scalar>::funPtrPhiPhi(P));
    phiPhiDx =
        (Basis::BasisHandler<
            typename LinearOperatorTraits<Derived>::Scalar>::funPtrPhiPhiDx(P));
    phiPhiDy =
        (Basis::BasisHandler<
            typename LinearOperatorTraits<Derived>::Scalar>::funPtrPhiPhiDy(P));
    phiTimesPhi =
        (Basis::BasisHandler<typename LinearOperatorTraits<Derived>::Scalar>::
             funPtrPhiTimesPhi(P));
    // vPhiScalVPhi = (Basis::BasisHandler<typename
    // LinearOperatorTraits<Derived>::Scalar>::funPtrVPhiScalVPhi(P));
    divPhiTimesDivPhi =
        (Basis::BasisHandler<typename LinearOperatorTraits<Derived>::Scalar>::
             funPtrDivPhiTimesDivPhi(P));

    curlPhiTimesCurlPhi =
                 (Basis::BasisHandler<typename LinearOperatorTraits<Derived>::Scalar>::
                      funPtrCurlPhiTimesCurlPhi(P));

    vecTimesCurlPhi =
                 (Basis::BasisHandler<typename LinearOperatorTraits<Derived>::Scalar>::
                                        funPtrVecTimesCurlPhi(P));
    vecTimesCurlPhi_single_basis =
                 (Basis::BasisHandler<typename LinearOperatorTraits<Derived>::Scalar>::
                                        funPtrVecTimesCurlPhi_single_basis(P));

    vecTimesCurlPhi_BSpline =
                 (Basis::BasisHandler<typename LinearOperatorTraits<Derived>::Scalar>::
                                        funPtrVecTimesCurlPhi_BSpline(P));

    mesh_ = std::make_shared<ClusterTree>();
    mesh_->init_ClusterTree(geom, M);
    mesh_->checkOrientation();
    return;
  };
  //////////////////////////////////////////////////////////////////////////////
  //    map2surface
  //////////////////////////////////////////////////////////////////////////////
  void map2surface(const ElementTreeNode& e, const Eigen::Vector2d& xi,
                   double w, SurfacePoint* surf_pt) const {
    Eigen::Vector2d st = e.llc_ + e.get_h() * xi;
    mesh_->get_geometry()[e.patch_].updateSurfacePoint(surf_pt, st, w, xi);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    Methods
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Compute all products of local shape functions on the unit square at
   * coordinates s,t, scale by w and add to intval.
   */
  void addScaledBasisInteraction(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s, const Eigen::Vector2d& t) const {
    phiTimesPhi(intval, w, s, t);
  }

  /**
   * \brief Compute all curl(phiphi(x)) curl(phiphi(y)) scalar products at the
   points x,y with tangential vectors x_f_dx, x_f_dy, y_f_dx and y_f_dy
   and normal vectors n_x, n_y.
   */
  void addScaledCurlBasisInteraction(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s, const Eigen::Vector2d& t,
      const Eigen::Vector3d x_f_dx, const Eigen::Vector3d x_f_dy,
      const Eigen::Vector3d y_f_dx, const Eigen::Vector3d y_f_dy) const {
    curlPhiTimesCurlPhi(intval, w, s, t,x_f_dx,x_f_dy,y_f_dx,y_f_dy);

  }

  /**
   * \brief Compute all vec times curl(phiphi(x)) scalar products at the
   point x with tangential vector x_f_dx, x_f_dy and normal vector n_x.
   */
  void addScaledVecTimesCurlBasisInteraction(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>* intval,
      const Eigen::Vector3d vec,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s,  const Eigen::Vector3d x_f_dx,
      const Eigen::Vector3d x_f_dy) const {
    vecTimesCurlPhi(intval, vec,w, s, x_f_dx, x_f_dy);

  }

  /**
   * \brief Compute one vec times curl(phiphi(x)) scalar product at the
   point x with tangential vector x_f_dx, x_f_dy and normal vector n_x.
   */
  void addScaledVecTimesCurlBasisInteraction_single_basis(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>* intval,
      const Eigen::Vector3d vec,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s,  const Eigen::Vector3d x_f_dx,
      const Eigen::Vector3d x_f_dy, int k) const {
    vecTimesCurlPhi_single_basis(intval, vec,w, s, x_f_dx, x_f_dy, k);

  }

  /**
   * \brief Compute one vec times curl(phiphi(x)) scalar product at the
   point x with tangential vector x_f_dx, x_f_dy and normal vector n_x.
   */
  void addScaledVecTimesCurlBasisInteraction_BSpline(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>* intval,
      const Eigen::Vector3d vec,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s,  const Eigen::Vector3d x_f_dx,
      const Eigen::Vector3d x_f_dy, double *weights, int *brnst_fcns,int num_fcns) const {
    vecTimesCurlPhi_BSpline(intval, vec,w, s, x_f_dx, x_f_dy, weights, brnst_fcns, num_fcns);

  }

  /**
   * \brief Compute all products of local shape functions on the unit square at
   * coordinates s,t.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                Eigen::Dynamic>
  basisInteraction(const Eigen::Vector2d& s, const Eigen::Vector2d& t) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, Eigen::Dynamic>
        intval(polynomial_degree_plus_one_squared,
               polynomial_degree_plus_one_squared);
    intval.setZero();
    phiTimesPhi(&intval, 1., s, t);
    return intval;
  }
  /**
   * \brief Compute all scalar products of vector valued local shape functions
   * on the surface points with reference coordinates s,t, scale by w and add to
   * intval.
   */
  void addScaledVectorBasisInteraction(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s, const Eigen::Vector2d& t,
      const Eigen::Vector3d x_f_dx, const Eigen::Vector3d x_f_dy,
      const Eigen::Vector3d y_f_dx, const Eigen::Vector3d y_f_dy) const {
    auto basis_interaction = basisInteraction(s, t);
    intval->block(0, 0, polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dx.dot(y_f_dx) * basis_interaction;
    intval->block(0, polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dx.dot(y_f_dy) * basis_interaction;
    intval->block(polynomial_degree_plus_one_squared, 0,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dy.dot(y_f_dx) * basis_interaction;
    intval->block(polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared,
                  polynomial_degree_plus_one_squared) +=
        w * x_f_dy.dot(y_f_dy) * basis_interaction;
  }
  /**
   * \brief Compute all scalar products of vector valued local shape functions
   * on the surface points with reference coordinates s,t.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                Eigen::Dynamic>
  vectorBasisInteraction(const Eigen::Vector2d& s, const Eigen::Vector2d& t,
                         const Eigen::Vector3d x_f_dx,
                         const Eigen::Vector3d x_f_dy,
                         const Eigen::Vector3d y_f_dx,
                         const Eigen::Vector3d y_f_dy) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, Eigen::Dynamic>
        intval(2 * polynomial_degree_plus_one_squared,
               2 * polynomial_degree_plus_one_squared);
    intval.setZero();
    addScaledVectorBasisInteraction(&intval, 1., s, t, x_f_dx, x_f_dy, y_f_dx,
                                    y_f_dy);
    return intval;
  }
  /**
   * \brief Compute all products of divergences of local shape functions on the
   * unit square at coordinates s,t, scale by w and add to intval.
   */
  void addScaledVectorBasisDivergenceInteraction(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s, const Eigen::Vector2d& t) const {
    divPhiTimesDivPhi(intval, w, s, t);
  }
  /**
   * \brief Compute all products of divergences of local shape functions on the
   * unit square at coordinates s,t.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                Eigen::Dynamic>
  vectorBasisDivergenceInteraction(const Eigen::Vector2d& s,
                                   const Eigen::Vector2d& t) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, Eigen::Dynamic>
        intval(2 * polynomial_degree_plus_one_squared,
               2 * polynomial_degree_plus_one_squared);
    intval.setZero();
    divPhiTimesDivPhi(&intval, 1., s, t);
    return intval;
  }
  /**
   * \brief Evaluate local shape functions on the unit square at coordinate s,
   * scale by w and add to intval.
   */
  void addScaledBasis(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, 1>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s) const {
    phiPhi(intval, w, s);
  }
  /**
   * \brief Evaluate local shape functions on the unit square at coordinate s.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                1>
  basis(const Eigen::Vector2d& s) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, 1>
        intval(polynomial_degree_plus_one_squared);
    intval.setZero();
    phiPhi(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate derivatives in x direction of local shape functions on the
   * unit square at coordinate s, scale by w and add to intval.
   */
  void addScaledBasisDx(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, 1>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s) const {
    phiPhiDx(intval, w, s);
  }
  /**
   * \brief Evaluate derivatives in x direction of local shape functions on the
   * unit square at coordinate s.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                1>
  basisDx(const Eigen::Vector2d& s) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, 1>
        intval(polynomial_degree_plus_one_squared);
    intval.setZero();
    phiPhiDx(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate derivatives in y direction of local shape functions on the
   * unit square at coordinate s, scale by w and add to intval.
   */
  void addScaledBasisDy(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, 1>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w,
      const Eigen::Vector2d& s) const {
    phiPhiDy(intval, w, s);
  }
  /**
   * \brief Evaluate derivatives in y direction of local shape functions on the
   * unit square at coordinate s.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                1>
  basisDy(const Eigen::Vector2d& s) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, 1>
        intval(polynomial_degree_plus_one_squared);
    intval.setZero();
    phiPhiDy(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate local shape functions on the unit interval at coordinate s,
   * scale by w and add to intval.
   */
  void addScaledBasis1D(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, 1>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w, double s) const {
    phi(intval, w, s);
  }
  /**
   * \brief Evaluate local shape functions on the unit interval at coordinate s.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                1>
  basis1D(double s) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, 1>
        intval(polynomial_degree + 1);
    intval.setZero();
    phi(&intval, 1., s);
    return intval;
  }
  /**
   * \brief Evaluate derivatives of local shape functions on the unit interval
   * at coordinate s, scale by w and add to intval.
   */
  void addScaledBasis1DDx(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, 1>* intval,
      typename LinearOperatorTraits<Derived>::Scalar w, double s) const {
    phiDx(intval, w, s);
  }
  /**
   * \brief Evaluate derivatives of local shape functions on the unit interval
   * at coordinate s.
   */
  Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
                1>
  basis1DDx(double s) const {
    Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                  Eigen::Dynamic, 1>
        intval(polynomial_degree + 1);
    intval.setZero();
    phiDx(&intval, 1., s);
    return intval;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  std::shared_ptr<ClusterTree> mesh_;
  Basis::funptr_phi<typename LinearOperatorTraits<Derived>::Scalar> phi;
  Basis::funptr_phidx<typename LinearOperatorTraits<Derived>::Scalar> phiDx;
  Basis::funptr_phiphi<typename LinearOperatorTraits<Derived>::Scalar> phiPhi;
  Basis::funptr_phiphidx<typename LinearOperatorTraits<Derived>::Scalar>
      phiPhiDx;
  Basis::funptr_phiphidy<typename LinearOperatorTraits<Derived>::Scalar>
      phiPhiDy;
  Basis::funptr_phitimesphi<typename LinearOperatorTraits<Derived>::Scalar>
      phiTimesPhi;
  // Basis::funptr_vphiscalvphi<typename LinearOperatorTraits<Derived>::Scalar>
  // vPhiScalVPhi;
  Basis::funptr_divphitimesdivphi<
      typename LinearOperatorTraits<Derived>::Scalar>
      divPhiTimesDivPhi;

  Basis::funptr_curlphitimescurlphi<
          typename LinearOperatorTraits<Derived>::Scalar>
          curlPhiTimesCurlPhi;

  Basis::funptr_vectimescurlphi<
                  typename LinearOperatorTraits<Derived>::Scalar>
                  vecTimesCurlPhi;

  Basis::funptr_vectimescurlphi_single_basis<
                  typename LinearOperatorTraits<Derived>::Scalar>
                  vecTimesCurlPhi_single_basis;

  Basis::funptr_vectimescurlphi_BSpline<
                  typename LinearOperatorTraits<Derived>::Scalar>
                  vecTimesCurlPhi_BSpline;

  int polynomial_degree;
  int polynomial_degree_plus_one_squared;
};  // namespace Bembel
}  // namespace Bembel
#endif
