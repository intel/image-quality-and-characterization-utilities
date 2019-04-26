#!/bin/bash
GIVEN_PARAMETERS=$@

function usage()
{
  echo "$(basename ${0}) [-t <tempdir>] [-l <logdir>] [-p <package-dir>] [-m <email-template>] [-b] [-s] [-k] [-f] [-g <GIT-project>] [-v <version>]"
  echo "   -t <tempdir>         overrides the directory where the build intermediates are stored"
  echo "   -l <logdir>          overrides the directory where the logfiles are stored"
  echo "   -p <packages-dir>    overrides the directory where the final build results are stored"
  echo "   -m <email-template>  points to a file that is used to collect build results, the contents of this file"
  echo "                        will be part of the mail send to the people also specified in this file"
  echo "                        To specify recipients, specify each recipient email-address with a line that starts"
  echo "                        with 'TO:' (whitespace seperated, only plain email-addresses)"
  echo "   -b                   only build with no test"
  echo "   -s                   'Smoketest-build' Builds the packages, but only runs a small smoketest for verification"
  echo "                        This is part of the 'build-on-commit' process executed by the CI-system"
  echo "   -v <version>         version number of this build"
  echo "   -k                   Klocwork build (no bsub, PAR or KERNEL_PAR)"
  echo "   -f                   'Fulltest-build' runs full test, useful for test taking many days to execute"
  echo "   -h                   this help"
  echo "   -i <patchid>         part of patch driver version, the patch changeid and patchset number"
  echo "   -d <build_type>      to identify whether this is a patchset build(un-merged) or a build from trunk(merged) (tip=false of true)"
  echo "   -n <build_number>    build_number in TeamCity"
  exit 1
}

function detect_os()
{
  MY_OS="$(uname -s)"
  case $MY_OS in
    Linux)
      #For now assume SLES11 environment
      os_env=Linux
      ;;
    CYGWIN*|MINGW*|MSYS*)
      os_env=Windows
      ;;
    *)
     echo "Error: Unsupported OS $MY_OS"
     os_env=Unsupported
     exit 1
     ;;
  esac

  echo "OS environment : $MY_OS"
  echo "Build category : $os_env"
}

function parse_args()
{
  temp_dir="";       temp_dir_opt=""
  log_dir="";        log_dir_opt=""
  pkg_dir="";        pkg_dir_opt=""
  mail_template="";  mail_template_opt=""
  version="";        version_opt=""
  git_project=""
  smoketest=""
  klocwork=""
  fulltest=""
  build_only=""
  regress_only=""
  patchid=""
  build_type=""
  build_number=0

  while getopts "t:l:p:m:v:g:i:d:n:skfbrh" OPTION; do
  case $OPTION in
      t)
    temp_dir="${OPTARG}"
    temp_dir_opt="-t"
    echo "temp_dir       : ${temp_dir}"
    ;;
      l)
    log_dir="${OPTARG}"
    log_dir_opt="-l"
    echo "log_dir        : ${log_dir}"
    ;;
      p)
    pkg_dir="${OPTARG}"
    pkg_dir_opt="-p"
    echo "pkg_dir        : ${pkg_dir}"
    ;;
      m)
    mail_template="${OPTARG}"
    mail_template_opt="-m"
    rm -f ${mail_template}
    echo "mail_template  : ${mail_template}"
    ;;
      v)
    version="${OPTARG}"
    version_opt="-v"
    echo "version        : ${version}"
    ;;
      g)
    git_project="${OPTARG}"
    echo "git_project    : ${git_project}"
    ;;
      i)
    patchid="${OPTARG}"
    echo "patchid        : ${patchid}"
    ;;
      d)
    build_type="${OPTARG}"
    if [ "${build_type}" == "true" ];then
        build_type="trunk"
    elif [ "${build_type}" == "false" ]; then
        build_type="patchset"
    fi
    echo "build_type     : ${build_type}"
    ;;
      n)
    build_number="${OPTARG}"
    echo "build_number   : ${build_number}"
    ;;
      b)
    build_only="-b"
    ;;
      r)
    regress_only="-r"
    ;;
      s)
    smoketest="-s"
    ;;
      k)
    klocwork="-k"
    ;;
      f)
    fulltest="-f"
    ;;
      h)
    usage
    ;;
      ?)
    echo "See help -h for usage"
    exit 1
    ;;
  esac
    done
    return 0
}

function clean_build_area()
{
  rm -rf   $build_dir
  mkdir -p $build_dir
}

function build_and_test_with_gcc()
{
  # Configure, Build and Test, but Exit immediately on error.
  cmake -DCMAKE_BUILD_TYPE=$1 $source_dir  &> $log_dir/cmake_$1.log  || exit $?
  cmake --build .  --config $1 -- -j        > $log_dir/build_$1.log  || exit $?
  ctest --verbose --build-config $1         > $log_dir/test_$1.log   || exit $?
}

function build_on_linux()
{
  if [ -e ../.slicenv/env/slicenv.sourceme.sh ]
  then
    echo "Calling slic start imaging"
    source ../.slicenv/env/slicenv.sourceme.sh
    slic start imaging
  else
    echo "Could not find slicenv.sourceme.sh script"
    exit 1
  fi

  clean_build_area
  pushd $build_dir

  set -x
  mkdir debug   && pushd debug   && build_and_test_with_gcc debug   && popd
  mkdir release && pushd release && build_and_test_with_gcc release && popd
  set +x

  error_code=0
}

function build_on_windows()
{
  clean_build_area
  pushd $build_dir

  # Configure, Build and Test, but Exit immediately on error.
  generator="Visual Studio 15 2017 Win64"
  build_cmd="cmake --build . --config"
  test_cmd="ctest --verbose --build-config"
  set -x
  cmake -G "$generator" $source_dir  &> $log_dir/cmake.log         || exit $?
  $build_cmd Debug   -- -maxcpucount  > $log_dir/build_debug.log   || exit $?
  $build_cmd Release -- -maxcpucount  > $log_dir/build_release.log || exit $?
  $test_cmd Debug                     > $log_dir/test_debug.log    || exit $?
  $test_cmd Release                   > $log_dir/test_release.log  || exit $?
  set +x

  error_code=0
}

function main()
{
  detect_os
  parse_args $GIVEN_PARAMETERS

  # Make sure that paths exists
  mkdir -p $temp_dir
  mkdir -p $log_dir

  # https://stackoverflow.com/questions/59895/get-the-source-directory-of-a-bash-script-from-within-the-script-itself
  source_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
  build_dir=$temp_dir/_build

  case ${os_env} in
    Linux|linux)
      build_on_linux
      ;;
    Windows|windows*)
      build_on_windows
      ;;
    *)
      echo "Error: Unsupported OS $os_env"
      error_code=9
      ;;
  esac

  popd
  echo "Error Code: ${error_code}"
  exit ${error_code}
}

main
