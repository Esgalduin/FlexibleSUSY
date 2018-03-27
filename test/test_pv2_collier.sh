#!/bin/sh

BASEDIR=$(dirname $0)
FSCONFIG="${BASEDIR}/../flexiblesusy-config"

src_FS="${BASEDIR}/test_pv2_run.cpp"
run_FS="${BASEDIR}/test_pv2_run.x"
src_CO="${BASEDIR}/test_pv2_run_collier.F90"
run_CO="${BASEDIR}/test_pv2_run_collier.x"
input="${BASEDIR}/test_pv2_collier_2point.in.dat"
out_FS="${BASEDIR}/test_pv2_collier_2point.FlexibleSUSY.out.dat"
out_CO="${BASEDIR}/test_pv2_collier_2point.COLLIER.out.dat"

[ -e "${run_FS}" -a "${run_FS}" -nt "${src_FS}" ] || {
    echo "${run_FS} must be rebuild, building it now ..."
    eval "$("${FSCONFIG}" --compile-cmd) ${src_FS} -o ${run_FS} $("${FSCONFIG}" --libs)"
    [ $? != 0 ] && {
        echo "Error: could not build ${run_FS}"
        exit 0
    }
}

# compile with
#
#  $ COLLIER_DIR=${HOME}/packages/COLLIER-1.2/
#
#  $ gfortran test/test_pv2_run_collier.F90 \
#      -o test/test_pv2_run_collier.x \
#      -J${COLLIER_DIR}/modules -L${COLLIER_DIR} -lcollier

[ -e "${run_CO}" -a "${run_CO}" -nt "${src_CO}" ] || {
    echo "${run_CO} must be rebuild, building it now ..."

    [ -z "${COLLIER_DIR}" ] && {
        echo "Please set COLLIER_DIR to the path of your COLLIER installation"
        exit 0
    }

    $("${FSCONFIG}" --fc) "${src_CO}" -o "${run_CO}" -J"${COLLIER_DIR}/modules" -L"${COLLIER_DIR}" -lcollier

    [ $? != 0 ] && {
        echo "Error: could not build ${run_FS}"
        exit 0
    }
}

command -v numdiff > /dev/null || {
    echo "numdiff is not installed"
    exit 0
}

"${run_FS}" "${input}" > "${out_FS}"
"${run_CO}" "${input}" > "${out_CO}"

# strip COLLIER header

sed -i -e '/^ *$/d' -e '/ *\*/d' "${out_CO}"

rel_error=4e-6
abs_error=1e-6

diff=$(numdiff --relative-tolerance="$rel_error" --absolute-tolerance="$abs_error" "${out_FS}" "${out_CO}")
diff=$(echo "$diff" | sed -e '/^ *#/d' | sed -e '/^+++/d')

if [ -n "$diff" ]; then
    echo "Error: difference between ${out_FS} and ${out_CO} larger that $rel_error"
    echo "$diff"
    echo ""
    echo "Test result: FAIL"
    exit_code=1
else
    echo ""
    echo "Test result: OK"
    exit_code=0
fi

exit $exit_code
