if [ ! -d "python-plots" ]; then
  mkdir python-plots
fi

if [ ! -d "PositivityFourMoments" ]; then
  mkdir PositivityFourMoments
fi
mv *positivity_shift_FourMoments_* PositivityFourMoments/

if [ ! -d "BiMaxwellianMoments" ]; then
  mkdir BiMaxwellianMoments
fi
mv *_BiMaxwellianMoments_* BiMaxwellianMoments/

if [ ! -d "Field" ]; then
  mkdir Field
fi
mv *-field_* Field/

if [ ! -d "Geometry" ]; then
  mkdir Geometry
fi

mv *-b_i* Geometry/
mv *-bcart* Geometry/
mv *-bmag* Geometry/
mv *-cmag* Geometry/
mv *-dxdz* Geometry/
mv *-dzdx* Geometry/
mv *-eps2* Geometry/
mv *-g_ij* Geometry/
mv *-gij* Geometry/
mv *-gxxj* Geometry/
mv *-gxyj* Geometry/
mv *-gxzj* Geometry/
mv *-gyyj* Geometry/
mv *-jacobgeo* Geometry/
mv *-jacobtot* Geometry/
mv *-mapc2p* Geometry/
mv *-mc2nu_pos.gkyl* Geometry/
mv *-nodes* Geometry/
mv *-normals* Geometry/

if [ ! -d "Slurmscripts" ]; then
  mkdir Slurmscripts
fi
mv *.out Slurmscripts/

if [ ! -d "Source" ]; then
  mkdir Source
fi
mv *_source_* Source/

if [ ! -d "Coll" ]; then
  mkdir Coll
fi
mv *_nu_sum_* Coll/
mv *_nu_prim_moms_* Coll/

if [ ! -d "PrimMoms" ]; then
  mkdir PrimMoms
fi
mv *_prim_moms_* PrimMoms/

if [ ! -d "M" ]; then
  mkdir M
fi
mv *_M* M/

if [ ! -d "Distributions" ]; then
  mkdir Distributions
fi
mv *-ion_[0-9]* Distributions/
mv *-elc_[0-9]* Distributions/

if [ ! -d "misc" ]; then
  mkdir misc
fi
mv *.gkyl misc/
mv *.json misc/