#ifndef fmssmn_lattice_constraints_hpp
#define fmssmn_lattice_constraints_hpp


#include "lattice_foreign_constraint.hpp"
#include "small_matrices.hpp"


namespace flexiblesusy {

#define decl_fmssmn_bc(name)						\
									\
extern "C" void name##_							\
(const Real& g1i, const Real& g2i, const Real& g3i,			\
 const Comp *Yui, const Comp *Ydi, const Comp *Yni, const Comp *Yei,	\
 const Real& m2Hui, const Real& m2Hdi,					\
 const Comp *m2Qi, const Comp *m2Ui, const Comp *m2Di,			\
 const Comp *m2Li, const Comp *m2Ni, const Comp *m2Ei,			\
 const Comp *Aui, const Comp *Adi, const Comp *Ani, const Comp *Aei,	\
 const Comp& M1i, const Comp& M2i, const Comp& M3i,			\
 const Real& vu, const Real& vd,					\
 const Real& scale0, const Real *x, const int& i,			\
 Real *row, Real *rhs);

decl_fmssmn_bc(fmssmn_mx)
decl_fmssmn_bc(fmssmn_higgs_masses)
decl_fmssmn_bc(fmssmn_gaugino_masses)
decl_fmssmn_bc(fmssmn_sfermion_masses)
decl_fmssmn_bc(fmssmn_trilinear_factors)
decl_fmssmn_bc(fmssmn_real_trilinear_factors)
decl_fmssmn_bc(fmssmn_ms)
decl_fmssmn_bc(fmssmn_gauge_couplings)
decl_fmssmn_bc(fmssmn_yude)
decl_fmssmn_bc(fmssmn_yn)
decl_fmssmn_bc(fmssmn_ewsb)


class Fmssmn_constraint_on_mx : public ForeignConstraint {
public:
    Fmssmn_constraint_on_mx() : ForeignConstraint(1) {}
    void operator()() {
	fmssmn_mx_(0,0,0,
		   nullptr,nullptr,nullptr,nullptr,
		   0,0,
		   nullptr,nullptr,nullptr,
		   nullptr,nullptr,nullptr,
		   nullptr,nullptr,nullptr,nullptr,
		   0,0,0,
		   0,0,
		   f->scl0, nullptr, 0,
		   rows.data(), &rhss(0));
	copy_rows();
    }
};

class Fmssmn_constraint_on_ms : public ForeignConstraint {
public:
    Fmssmn_constraint_on_ms() : ForeignConstraint(1) {}
    void operator()() {
	fmssmn_ms_(0,0,0,
		   nullptr,nullptr,nullptr,nullptr,
		   0,0,
		   nullptr,nullptr,nullptr,
		   nullptr,nullptr,nullptr,
		   nullptr,nullptr,nullptr,nullptr,
		   0,0,0,
		   vu,vd,
		   f->scl0, &x()[0], 0,
		   rows.data(), &rhss(0));
	copy_rows();
    }
    Real vu, vd;
};

class Fmssmn_constraint_on_gauge_couplings : public ForeignConstraint {
public:
    Fmssmn_constraint_on_gauge_couplings() : ForeignConstraint(3) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 3; i++) {
	    fmssmn_gauge_couplings_(g1,g2,g3,
				    nullptr,nullptr,nullptr,nullptr,
				    0,0,
				    nullptr,nullptr,nullptr,
				    nullptr,nullptr,nullptr,
				    nullptr,nullptr,nullptr,nullptr,
				    0,0,0,
				    0,0,
				    f->scl0, nullptr, i,
				    rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    Real g1, g2, g3;
};

class Fmssmn_constraint_on_yude : public ForeignConstraint {
public:
    Fmssmn_constraint_on_yude() : ForeignConstraint(54) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 54; i++) {
	    fmssmn_yude_(0,0,0,
			 Yu.data(),Yd.data(),nullptr,Ye.data(),
			 0,0,
			 nullptr,nullptr,nullptr,
			 nullptr,nullptr,nullptr,
			 nullptr,nullptr,nullptr,nullptr,
			 0,0,0,
			 0,0,
			 f->scl0, nullptr, i,
			 rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    CM33 Yu, Yd, Ye;
};


class Fmssmn_constraint_on_ewsb : public ForeignConstraint {
public:
    Fmssmn_constraint_on_ewsb() : ForeignConstraint(4) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 4; i++) {
	    fmssmn_ewsb_(0,0,0,
			 nullptr,nullptr,nullptr,nullptr,
			 0,0,
			 nullptr,nullptr,nullptr,
			 nullptr,nullptr,nullptr,
			 nullptr,nullptr,nullptr,nullptr,
			 0,0,0,
			 vu,vd,
			 f->scl0, &x()[0], i,
			 rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    Real vu, vd;
};

class Fmssmn_constraint_on_higgs_masses : public ForeignConstraint {
public:
    Fmssmn_constraint_on_higgs_masses() : ForeignConstraint(2) {}
    void operator()() {
	Real m2Hu = (1-f->a)*m2Hu_ini + f->a*m2Hu_fin;
	Real m2Hd = (1-f->a)*m2Hd_ini + f->a*m2Hd_fin;
	rows.transposeInPlace();
	for (size_t i = 0; i < 2; i++) {
	    fmssmn_higgs_masses_(0,0,0,
				 nullptr,nullptr,nullptr,nullptr,
				 m2Hu,m2Hd,
				 nullptr,nullptr,nullptr,
				 nullptr,nullptr,nullptr,
				 nullptr,nullptr,nullptr,nullptr,
				 0,0,0,
				 0,0,
				 f->scl0, nullptr, i,
				 rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    Real m2Hu_ini, m2Hu_fin;
    Real m2Hd_ini, m2Hd_fin;
};

class Fmssmn_constraint_on_gaugino_masses : public ForeignConstraint {
public:
    Fmssmn_constraint_on_gaugino_masses() : ForeignConstraint(6) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 6; i++) {
	    fmssmn_gaugino_masses_(0,0,0,
				   nullptr,nullptr,nullptr,nullptr,
				   0,0,
				   nullptr,nullptr,nullptr,
				   nullptr,nullptr,nullptr,
				   nullptr,nullptr,nullptr,nullptr,
				   M1,M2,M3,
				   0,0,
				   f->scl0, nullptr, i,
				   rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    Comp M1, M2, M3;
};

class Fmssmn_constraint_on_sfermion_masses : public ForeignConstraint {
public:
    Fmssmn_constraint_on_sfermion_masses() : ForeignConstraint(54) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 54; i++) {
	    fmssmn_sfermion_masses_(0,0,0,
				    nullptr,nullptr,nullptr,nullptr,
				    0,0,
				    m2Q.data(),m2U.data(),m2D.data(),
				    m2L.data(),m2N.data(),m2E.data(),
				    nullptr,nullptr,nullptr,nullptr,
				    0,0,0,
				    0,0,
				    f->scl0, nullptr, i,
				    rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    CM33 m2Q, m2U, m2D, m2L, m2N, m2E;
};

class Fmssmn_constraint_trilinear_factors : public ForeignConstraint {
public:
    Fmssmn_constraint_trilinear_factors() : ForeignConstraint(72) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 72; i++) {
	    fmssmn_trilinear_factors_(0,0,0,
				      nullptr,nullptr,nullptr,nullptr,
				      0,0,
				      nullptr,nullptr,nullptr,
				      nullptr,nullptr,nullptr,
				      Au.data(),Ad.data(),An.data(),Ae.data(),
				      0,0,0,
				      0,0,
				      f->scl0, nullptr, i,
				      rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    CM33 Au, Ad, An, Ae;
};

class Fmssmn_constraint_real_trilinear_factors : public ForeignConstraint {
public:
    Fmssmn_constraint_real_trilinear_factors() : ForeignConstraint(72) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 72; i++) {
	    fmssmn_real_trilinear_factors_(0,0,0,
				      nullptr,nullptr,nullptr,nullptr,
				      0,0,
				      nullptr,nullptr,nullptr,
				      nullptr,nullptr,nullptr,
				      Au.data(),Ad.data(),An.data(),Ae.data(),
				      0,0,0,
				      0,0,
				      f->scl0, nullptr, i,
				      rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    CM33 Au, Ad, An, Ae;	// imaginary parts are unused
};

class Fmssmn_constraint_on_yn : public ForeignConstraint {
public:
    Fmssmn_constraint_on_yn() : ForeignConstraint(18) {}
    void operator()() {
	rows.transposeInPlace();
	for (size_t i = 0; i < 18; i++) {
	    fmssmn_yn_(0,0,0,
		       nullptr,nullptr,Yn.data(),nullptr,
		       0,0,
		       nullptr,nullptr,nullptr,
		       nullptr,nullptr,nullptr,
		       nullptr,nullptr,nullptr,nullptr,
		       0,0,0,
		       0,0,
		       f->scl0, nullptr, i,
		       rows.col(i).data(), &rhss(i));
	}
	rows.transposeInPlace();
	copy_rows();
    }
    CM33 Yn;
};

}

#endif // fmssmn_lattice_constraints_hpp
