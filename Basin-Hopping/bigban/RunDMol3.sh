#!/bin/sh 
    # ---------------------------------------------------------------------------
    #
    # Script for stand-alone DMol3 execution.
    #
    # This program and all subroutines, data, and files used by it
    # are protected by Copyright and hence may not be used, copied,
    # modified, transmitted, inspected, or executed by any means including
    # the use of electronic data processing equipment, xerography, or
    # any other methods without the express written permission of the
    # copyright holder.
    #
    # Copyright (C) 2004 Accelrys, Inc. All rights reserved
    #
    # ***************************************************************************
    
    usage () {
    SCRIPT_NAME=`basename $0`
    echo ""
    echo "DESCRIPTION"
    echo "  $SCRIPT_NAME script runs DMol3 as a standalone program"
    echo "  in serial or parallel mode."
    echo ""
    echo "  Typical usage:  "
    echo "     $SCRIPT_NAME [-h] [-nodelete] [-np <number of processors>] <basename>"
    echo ""
    echo "RUN OPTIONS"
    echo "     -h        - displays this help text. "
    echo "     -nodelete - retain all job and scratch files. "
    echo "     -np       - run in parallel on \"number of processors\". "
    echo "                 Serial version of DMol3 is executed when the "
    echo "                 -np option is missing. If the -np option is 1, then "
    echo "                 parallel DMol3 is executed on one CPU."
    echo ""
    echo "ARGUMENTS "
    echo "     basename  - name of the DMol3 job. Input is read from \"basename\".input "
    echo "                 and output is written to \"basename\".outmol"
    echo ""
    echo "  Note that \"basename\" should consist of alphanumeric characters"
    echo "  and can not contain any spaces."
    echo ""
    }
    
    tru64lim() {
    #increase user limits if this tru64 o/s
    os=`uname`
    if [ "$os" = "OSF1" ]; then
    # set stack size
      slimit=`/sbin/sysconfig -q proc  max-per-proc-stack-size | grep -v proc: | awk ' { print $3/1024 } '`
      ulimit -s $slimit
    # set data size
      dlimit=`/sbin/sysconfig -q proc  max-per-proc-data-size | grep -v proc: | awk ' { print $3/1024 } '`
      ulimit -d $dlimit
    fi
    }
    
    dmol3_parallel() {
      echo "executing parallel DMol3 on $nproc processors"
      os=`uname`
      if [ "$os" = "IRIX" ] || [ "$os" = "IRIX64" ]; then
         MPI_DSM_PPM=1
         export MPI_DSM_PPM
         MPC_GANG=OFF
         export MPC_GANG
         mpirun -np $nproc $server $rootname
      elif [ "$os" = "Linux" ]; then
          if [ $arc = "cluster" ]; then
              echo " "
              echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
              echo " HP-MPI runs DMol3 on a cluster machine on $nproc processors "
              echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
   
              envVars="-env DMOL3_DATA -env DMOL_TMP"
              $MS_INSTALL_ROOT/share/bin/create_appfile.sh -np $nproc $envVars -command "$server $rootname" -list $allnodes
   
              $MPI_COMMAND 
          else
              echo " "
              echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
              echo " HP-MPI runs DMol3 on an non-cluster SMP machine on $nproc processors "
              echo "    command: $MPI_COMMAND  "
              echo " processors: $nproc      "
              echo " executable: $server     "
              echo "   seedname: $rootname   "
              echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
              echo " "
              $MPI_COMMAND  $nproc $server "$rootname" $DMOL3_DATA < /dev/null
          fi
          
      elif [ "$os" = "OSF1" ]; then
        dmpirun -np $nproc $server $rootname
      else
        echo "Parallel runs are not supported on operating system $os by this script"
        exit 1
    
      fi
    }
    
    cluster_or_SMP() {
    
    arc=SMP
    echo "$MPI_COMMAND " | grep -e "-f"; retcode=$?
    if [ $retcode = 0 ]; then
        arc=cluster
    fi
    }
    
    # This is the "standard" name for the DMOL_TMP directory if the
    # user has chosen to have it "shared" within the job directory.
    # (This variable is read-only; its value is derived from the
    # installation scripts, but there is no mechanism for getting
    # the value directly.)
    DMOL_TMP_IN_JOB_FOLDER="./tmp"
 
    # "True" if we should attempt to clean up the tmpdir at
    # the end of the job.
    DMOL_TMP_IS_REMOVABLE=
 
    # Set the TMPDIR environment variable.
    # Linux cluster sometimes requires this to be in the job folder so
    # that it is shared between the nodes; at other times, we can simply
    # use the node-specific standard tmp.
    # NB This uses the variable $arc, so cluster_or_smp() must be called
    # before this.
    setDMOL_TMP() {
       if [ "$GATEWAY_TMP" != "" ]; then
          DMOL_TMP=$GATEWAY_TMP
          # Bug 04356vymq01: if $tmpdir is a relative, the DMol executable can't use it.
          # Rather than fix the executable, we ensure that the tmp dir is an absolute path:
          # we identify relative paths as not starting with a "/" and simply add the
          # job directory to the start.
          firstChar=`expr substr $DMOL_TMP 1 1`
          if [ "$firstChar" != "/" ]; then
              DMOL_TMP=`pwd`/$DMOL_TMP
              # Mark the tmp dir as available for deletion.
              DMOL_TMP_IS_REMOVABLE="true"
          fi
       else
          DMOL_TMP=$DMOL_TMP_IN_JOB_FOLDER
          DMOL_TMP_IS_REMOVABLE="true"
       fi
 
       # If tmpdir already exists, then we should leave it alone
       if [ -d $DMOL_TMP ]; then
          DMOL_TMP_IS_REMOVABLE=
       fi
    }
  
    clean_scratch() {
       #Remove any scratch files
       echo "cleaning scratch files"
       for ext in rot tpotl tmesh tpsmx optabs basis \
            grad inatm incoor opt outatom pchk \
               symdec symdv sym torder tpdiis tpdiisk tplev1 \
               tplev2 fwv amat prf monitor tphmx tpeig
       do
          [ -f $rootname.$ext ] && rm $rootname.$ext
       done
    
       rm -f *.pid
    }
    
    nproc=0
    nodelete=0
    while [ $# -gt 0 ]
    do
        case $1 in
            -np)      shift; nproc=$1;;
            -nodelete) nodelete=1;;	   
            -h)       usage; exit 0 ;;
            *)        rootname=$1;;
        esac
        shift
    done
    
    if [ "$rootname" = "" ]; then
      usage
      exit 1
    fi
    
    # Location of Materials Studio installation
    MS_INSTALL_ROOT=/home/MSI/MSI4.0
    export MS_INSTALL_ROOT
    # Location of the machines.LINUX file (cluster regime)
    allnodes=$MS_INSTALL_ROOT/share/data/machines.LINUX
    export allnodes
    
    #server executable
    if [ "$nproc" = "0" ]; then
      server="$MS_INSTALL_ROOT/DMol3/bin/dmol3.exe"
    else
      server="$MS_INSTALL_ROOT/DMol3/bin/dmol3_mpi.exe"
    fi
    
    #setup licensing
    if [ -f $MS_INSTALL_ROOT/share/license/data/lic_setup.sh ]; then
             eval `$MS_INSTALL_ROOT/share/license/data/lic_setup.sh $MS_INSTALL_ROOT/ -s sh`
    fi
    
    #setup ms execution environment
    if [ -f $MS_INSTALL_ROOT/share/bin/ms_setup.sh ]; then
             eval `$MS_INSTALL_ROOT/share/bin/ms_setup.sh $MS_INSTALL_ROOT/ -s sh`
    fi
    
    DMOL3_DATA=$MS_INSTALL_ROOT/Data/Resources/Quantum/DMol3
    export DMOL3_DATA
    
    #if this is tru64 o/s, then increase limits to the max allowed
    tru64lim
    
    setDMOL_TMP
    export DMOL_TMP
    if [ "$DMOL_TMP_IS_REMOVABLE" = "true" -a ! -d $DMOL_TMP ]; then
      mkdir -p $DMOL_TMP
    fi
    
    if [ "$nproc" = "0" ]; then
      #execute serial version
      echo "Executing serial DMol3"
      echo "$server $rootname"
      $server $rootname
    else 
      #execute parallel run
      cluster_or_SMP
      dmol3_parallel
    fi
 
    if [ "$DMOL_TMP_IS_REMOVABLE" = "true" ]; then
      echo "removing scratch directories '$DMOL_TMP'"
      /bin/rm -rf $DMOL_TMP
    fi
    
    if [ "$nodelete" = "0" ]; then
        clean_scratch
    fi  
    echo "DMol3 has completed"
    
    exit 0
