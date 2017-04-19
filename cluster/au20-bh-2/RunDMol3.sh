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
   # Copyright (C) 2007 Accelrys, Inc. All rights reserved
   #
   # ***************************************************************************
   MS_INSTALL_ROOT=/home/liyd/opt/Accelrys/MaterialsStudio7.0
   export MS_INSTALL_ROOT
   server=DMol3
   $MS_INSTALL_ROOT/share/bin/runMSserver.sh $server "$@"
   if [ $? != 6 ]; then
       exit $?
   fi

   cat $MS_INSTALL_ROOT/etc/DMol3/bin/RunDMol3.Readme
