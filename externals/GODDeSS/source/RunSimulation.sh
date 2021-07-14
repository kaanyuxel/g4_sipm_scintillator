#!/bin/bash


#/*
# * author:      Erik Dietz-Laursonn
# * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
# * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
# */


#--- program options to be specified: ---#

### Path to init-file containing program options to be passed to the program (it can be an existing init-file or the path where to save the used program options)
### Default: "$OutputDirectory/RunSimulation.init"
# InitFilePath=""
# InitFileName=""

### Path to macro-file containing command to be passed to the program
# MacroFilePath=""
# MacroFileName=""

### Number of events to be simulated in batch mode. If NumberOfEvents is set and greater than 0, the program will be running in batch mode (without visualisation, i.e. faster).
### Default: 0
# NumberOfEvents=1

### Use a list with photons to define the primary particles instead of using the normal particle source. This list has to be specified with "ParticleSourceMacroPath".
### Default: false
# UsePhotonList=true

### Path to input file for the particle source. This can either be a macro with commands for the General Particle Source (for UsePhotonList=false) or a gziped (filename.gz) list with photons (for UsePhotonList=true).
ParticleSourceMacroPath="$SIMDIR/macros"
ParticleSourceMacroName="GeneralParticleSource.mac"

### Use different dimensions (in mm) for the scintillator tile than the hardcoded ones
TileDimensions="(350,5,350)"

### Use the dimensions specified above to define the region in which the primary particles will be created (only if UsePhotonList=false)
### Default: false
UseTileDimensionsForParticleSource=true

### Search for overlaps in the setup
### Default: false
# SearchOverlaps=true

### Path to the directory containing the material property files
MaterialPropertiesPath="$SIMDIR/MaterialProperties"

### Path to material property file to be used for the scintillator material
ScintillatorMaterialPropertiesPath="$MaterialPropertiesPath/Scintillator"
# ScintillatorMaterialPropertiesName="Saint-Gobain_BC-404.properties"
ScintillatorMaterialPropertiesName="Saint-Gobain_BC-408.properties"
# ScintillatorMaterialPropertiesName="Saint-Gobain_BC-452_2perCentLead.properties"
# ScintillatorMaterialPropertiesName="Saint-Gobain_BC-452_5perCentLead.properties"
# ScintillatorMaterialPropertiesName="Saint-Gobain_BC-452_10perCentLead.properties"

### Path to material property file to be used for the twin tile reflector material
TwinTileReflectorMaterialPropertiesPath="$MaterialPropertiesPath/Scintillator"
TwinTileReflectorMaterialPropertiesName="Wrapping_Aluminum.properties"

### Path to material property file to be used for the wrapping material
WrappingMaterialPropertiesPath="$MaterialPropertiesPath/Scintillator"
# WrappingMaterialPropertiesName="Wrapping_Aluminum.properties"
WrappingMaterialPropertiesName="Wrapping_Teflon.properties"
# WrappingMaterialPropertiesName="Wrapping_Saint-Gobain_BC-620.properties"

### Path to material property file to be used for the material of light-guiding fibres
LGFibreMaterialPropertiesPath="$MaterialPropertiesPath/Fibre"
LGFibreMaterialPropertiesName="Saint-Gobain_BCF-98_round_1mm.properties"
# LGFibreMaterialPropertiesName="Saint-Gobain_BCF-98_round_2mm.properties"
# LGFibreMaterialPropertiesName="Saint-Gobain_BCF-98_quadratic_1mm.properties"
# LGFibreMaterialPropertiesName="Saint-Gobain_BCF-98_quadratic_2mm.properties"

### Path to material property file to be used for the material of WLS fibres
WLSFibreMaterialPropertiesPath="$MaterialPropertiesPath/Fibre"
WLSFibreMaterialPropertiesName="Saint-Gobain_BCF-92_round_1mm.properties"
# WLSFibreMaterialPropertiesName="Saint-Gobain_BCF-92_round_2mm.properties"
# WLSFibreMaterialPropertiesName="Saint-Gobain_BCF-92_quadratic_1mm.properties"
# WLSFibreMaterialPropertiesName="Saint-Gobain_BCF-92_quadratic_2mm.properties"

### Path to material property file to be used for the material of scintillating fibres
# ScintiFibreMaterialPropertiesPath="$MaterialPropertiesPath/Fibre"
# ScintiFibreMaterialPropertiesName=""

### Path to material property file to be used for the optical cement material
OpticalCementMaterialPropertiesPath="$MaterialPropertiesPath/OpticalCement"
OpticalCementMaterialPropertiesName="Saint-Gobain_BC-600.properties"
# OpticalCementMaterialPropertiesName="air.properties"

### Path to the output directory. (This directory must exisit!)
### Default: the directory in which the program runs
# OutputDirectory=""

### Name of the output files of the simulation.
### Default: "outputSimulation.data"
# SimulationOutputFileName="outputSimulation.data"

### Phrase to be added to the output file names
# AddStringToOutputFiles=""

### Verbosity
### Default: 0
# ControlVerbosity=1
# RunVerbosity=1
# EventVerbosity=1
# TrackingVerbosity=1

### Should optical photons be tracked  when simulating events?
### Default: true
# TrackOpticalPhotons=false



#--- process the program options ---#



if [ "$GODDESS" == "" ]; then
  echo
  echo "The environment variable \"GODDESS\" has not been specified!"
  echo

  exit 1
fi

if [ "$BUILDDIR" == "" ]; then
  echo
  echo "The environment variable \"BUILDDIR\" has not been specified!"
  echo

  exit 1
fi

if [ "$SIMDIR" == "" ]; then
  echo
  echo "The environment variable \"SIMDIR\" has not been specified!"
  echo

  exit 1
fi



# process the options
OtherOptionsGiven=false



if [ ! "$InitFilePath/$InitFileName" == "/" ] ; then
  if [ -f $InitFilePath/$InitFileName ] ; then
    InitFileExists=true
  else
    # abort the script
    echo "The specified init-file is not existing. The script is aborted"
    exit
  fi
else
  InitFileExists=false
fi


if [ ! "$OutputDirectory" == "" ] ; then
  if $InitFileExists ; then
    # abort the script
    echo "An existing init-file as well as other options are specified. The script is aborted"
    exit
  fi

  OutputDirectoryString="--outDir $OutputDirectory"
else
  if $InitFileExists && [ ! "`grep "\-\-outDir" $InitFilePath/$InitFileName`" == "" ] ; then
    OutputDirectory=`grep "\-\-outDir" $InitFilePath/$InitFileName`
    OutputDirectory=`echo $OutputDirectory | awk -F" " '{print $2}'`
  else
    OutputDirectory="$BUILDDIR"
  fi

  OutputDirectoryString=""
fi

if [ -e $OutputDirectory/_inputFiles/ ] ; then
  rm -rf $OutputDirectory/_inputFiles
fi
mkdir -p $OutputDirectory/_inputFiles/



if $InitFileExists ; then
  cp $InitFilePath/$InitFileName $OutputDirectory/_inputFiles/
else
  InitFilePath="$OutputDirectory"
  InitFileName="RunSimulation.init"

  if [ -f $InitFilePath/$InitFileName ] ; then
    rm $InitFilePath/$InitFileName
  fi
fi
InitFile="$OutputDirectory/_inputFiles/$InitFileName"
InitFileString="--init $OutputDirectory/_inputFiles/$InitFileName"


if [ ! "$MacroFilePath/$MacroFileName" == "/" ] ; then
  if [ ! -f $MacroFilePath/$MacroFileName ] ; then
    echo "\"$MacroFilePath/$MacroFileName\" does not exist!"
    exit 1
  fi

  cp $MacroFilePath/$MacroFileName $OutputDirectory/_inputFiles/
  MacroFileString="--macro $OutputDirectory/_inputFiles/$MacroFileName"

  OtherOptionsGiven=true
else
  MacroFileString=""
fi


if [ ! "$NumberOfEvents" == "" ] && [ $NumberOfEvents -gt 0 ] ; then
  RunBatchModeString="--batch $NumberOfEvents"

  OtherOptionsGiven=true
else
  RunBatchModeString=""
fi


if [ ! "$UsePhotonList" == "" ] && $UsePhotonList ; then
  UsePhotonListString="--usePhotonList"

  OtherOptionsGiven=true
else
  UsePhotonListString=""
fi


if [ ! "$ParticleSourceMacroPath/$ParticleSourceMacroName" == "/" ] ; then
  if [ ! -f $ParticleSourceMacroPath/$ParticleSourceMacroName ] ; then
    echo "\"$ParticleSourceMacroPath/$ParticleSourceMacroName\" does not exist!"
    exit 1
  fi

  cp $ParticleSourceMacroPath/$ParticleSourceMacroName $OutputDirectory/_inputFiles/
  ParticleSourceMacroPathString="--particleSourceInput $OutputDirectory/_inputFiles/$ParticleSourceMacroName"

  OtherOptionsGiven=true
else
  ParticleSourceMacroPathString=""
fi


if [ ! "$TileDimensions" == "" ] ; then
  TileDimensionsString="--tileDims $TileDimensions"

  OtherOptionsGiven=true
else
  TileDimensionsString=""
fi


if [ ! "$UseTileDimensionsForParticleSource" == "" ] && [ ! "$TileDimensions" == "" ] && $UseTileDimensionsForParticleSource && ( [ "$UsePhotonList" == "" ] || ! $UsePhotonList ) ; then
  ShapeSearchString=".*/gps/plane/shape.*"
  NewShapeLine="\t/gps/plane/shape\t\trect\t\t# Set the shape of the plane in which the primary particles will be created. Possible values are \"circle\", \"rect\" or \"hexagon\". (default: circle)"
  sed "s!$ShapeSearchString!$NewShapeLine!" -i "$OutputDirectory/_inputFiles/$ParticleSourceMacroName"

  XEdgeLength=`echo $TileDimensions | awk -F"," '{print $1}'`
  XEdgeLength=`echo $XEdgeLength | awk -F"(" '{print $2}'`

  XEdgeLengthSearchString=".*/gps/plane/a.*"
  NewXEdgeLengthLine="\t/gps/plane/a\t\t\t$XEdgeLength mm\t\t# If the shape of the plane in which the primary particles will be created is \"rect\", set its edge length in x direction. (default: 1 mm)"
  sed "s!$XEdgeLengthSearchString!$NewXEdgeLengthLine!" -i "$OutputDirectory/_inputFiles/$ParticleSourceMacroName"

  ZEdgeLength=`echo $TileDimensions | awk -F"," '{print $3}'`
  ZEdgeLength=`echo $ZEdgeLength | awk -F")" '{print $1}'`

  ZEdgeLengthSearchString=".*/gps/plane/b.*"
  NewZEdgeLengthLine="\t/gps/plane/b\t\t\t$ZEdgeLength mm\t\t# If the shape of the plane in which the primary particles will be created is \"rect\", set its edge length in y direction. (default: 1 mm)"
  sed "s!$ZEdgeLengthSearchString!$NewZEdgeLengthLine!" -i "$OutputDirectory/_inputFiles/$ParticleSourceMacroName"
fi


if [ ! "$SearchOverlaps" == "" ] && $SearchOverlaps ; then
  SearchOverlapsString="--overlap"

  OtherOptionsGiven=true
else
  SearchOverlapsString=""
fi


if [ ! "$ScintillatorMaterialPropertiesName" == "" ] ; then
  if [ ! -f $ScintillatorMaterialPropertiesPath/$ScintillatorMaterialPropertiesName ] ; then
    echo "\"$ScintillatorMaterialPropertiesPath/$ScintillatorMaterialPropertiesName\" does not exist!"
    exit 1
  fi

  cp $ScintillatorMaterialPropertiesPath/$ScintillatorMaterialPropertiesName $OutputDirectory/_inputFiles/
  ScintillatorMaterialPropertiesString="--scinti $OutputDirectory/_inputFiles/$ScintillatorMaterialPropertiesName"

  OtherOptionsGiven=true
else
  ScintillatorMaterialPropertiesString=""
fi


if [ ! "$TwinTileReflectorMaterialPropertiesName" == "" ] ; then
  if [ ! -f $TwinTileReflectorMaterialPropertiesPath/$TwinTileReflectorMaterialPropertiesName ] ; then
    echo "\"$TwinTileReflectorMaterialPropertiesPath/$TwinTileReflectorMaterialPropertiesName\" does not exist!"
    exit 1
  fi

  cp $TwinTileReflectorMaterialPropertiesPath/$TwinTileReflectorMaterialPropertiesName $OutputDirectory/_inputFiles/
  TwinTileReflectorMaterialPropertiesString="--twin $OutputDirectory/_inputFiles/$TwinTileReflectorMaterialPropertiesName"

  OtherOptionsGiven=true
else
  TwinTileReflectorMaterialPropertiesString=""
fi


if [ ! "$WrappingMaterialPropertiesName" == "" ] ; then
  if [ ! -f $WrappingMaterialPropertiesPath/$WrappingMaterialPropertiesName ] ; then
    echo "\"$WrappingMaterialPropertiesPath/$WrappingMaterialPropertiesName\" does not exist!"
    exit 1
  fi

  cp $WrappingMaterialPropertiesPath/$WrappingMaterialPropertiesName $OutputDirectory/_inputFiles/
  WrappingMaterialPropertiesString="--wrapping $OutputDirectory/_inputFiles/$WrappingMaterialPropertiesName"

  OtherOptionsGiven=true
else
  WrappingMaterialPropertiesString=""
fi


if [ ! "$LGFibreMaterialPropertiesName" == "" ] ; then
  if [ ! -f $LGFibreMaterialPropertiesPath/$LGFibreMaterialPropertiesName ] ; then
    echo "\"$LGFibreMaterialPropertiesPath/$LGFibreMaterialPropertiesName\" does not exist!"
    exit 1
  fi

  cp $LGFibreMaterialPropertiesPath/$LGFibreMaterialPropertiesName $OutputDirectory/_inputFiles/
  LGFibreMaterialPropertiesString="--lgFibre $OutputDirectory/_inputFiles/$LGFibreMaterialPropertiesName"

  OtherOptionsGiven=true
else
  LGFibreMaterialPropertiesString=""
fi


if [ ! "$WLSFibreMaterialPropertiesName" == "" ] ; then
  if [ ! -f $WLSFibreMaterialPropertiesPath/$WLSFibreMaterialPropertiesName ] ; then
    echo "\"$WLSFibreMaterialPropertiesPath/$WLSFibreMaterialPropertiesName\" does not exist!"
    exit 1
  fi

  cp $WLSFibreMaterialPropertiesPath/$WLSFibreMaterialPropertiesName $OutputDirectory/_inputFiles/
  WLSFibreMaterialPropertiesString="--wlsFibre $OutputDirectory/_inputFiles/$WLSFibreMaterialPropertiesName"

  OtherOptionsGiven=true
else
  WLSFibreMaterialPropertiesString=""
fi


if [ ! "$ScintiFibreMaterialPropertiesName" == "" ] ; then
  if [ ! -f $ScintiFibreMaterialPropertiesPath/$ScintiFibreMaterialPropertiesName ] ; then
    echo "\"$ScintiFibreMaterialPropertiesPath/$ScintiFibreMaterialPropertiesName\" does not exist!"
    exit 1
  fi

  cp $ScintiFibreMaterialPropertiesPath/$ScintiFibreMaterialPropertiesName $OutputDirectory/_inputFiles/
  ScintiFibreMaterialPropertiesString="--scintiFibre $OutputDirectory/_inputFiles/$ScintiFibreMaterialPropertiesName"

  OtherOptionsGiven=true
else
  ScintiFibreMaterialPropertiesString=""
fi


if [ ! "$OpticalCementMaterialPropertiesName" == "" ] ; then
  if [ ! -f $OpticalCementMaterialPropertiesPath/$OpticalCementMaterialPropertiesName ] ; then
    echo "\"$OpticalCementMaterialPropertiesPath/$OpticalCementMaterialPropertiesName\" does not exist!"
    exit 1
  fi

  cp $OpticalCementMaterialPropertiesPath/$OpticalCementMaterialPropertiesName $OutputDirectory/_inputFiles/
  OpticalCementMaterialPropertiesString="--cement $OutputDirectory/_inputFiles/$OpticalCementMaterialPropertiesName"

  OtherOptionsGiven=true
else
  OpticalCementMaterialPropertiesString=""
fi


if [ ! "$SimulationOutputFileName" == "" ] ; then
  SimulationOutputFileNameString="--outFile $SimulationOutputFileName"

  OtherOptionsGiven=true
else
  SimulationOutputFileNameString=""
fi


if [ ! "$AddStringToOutputFiles" == "" ] ; then
  AddStringToOutputFilesString="--add $AddStringToOutputFiles"

  OtherOptionsGiven=true
else
  AddStringToOutputFilesString=""
fi


if [ ! "$ControlVerbosity" == "" ] && [ $ControlVerbosity -gt 0 ] ; then
  ControlVerbosityString="--controlVerbose $ControlVerbosity"

  OtherOptionsGiven=true
else
  ControlVerbosityString=""
fi


if [ ! "$RunVerbosity" == "" ] && [ $RunVerbosity -gt 0 ] ; then
  RunVerbosityString="--runVerbose $RunVerbosity"

  OtherOptionsGiven=true
else
  RunVerbosityString=""
fi


if [ ! "$EventVerbosity" == "" ] && [ $EventVerbosity -gt 0 ] ; then
  EventVerbosityString="--eventVerbose $EventVerbosity"

  OtherOptionsGiven=true
else
  EventVerbosityString=""
fi


if [ ! "$TrackingVerbosity" == "" ] && [ $TrackingVerbosity -gt 0 ] ; then
  TrackingVerbosityString="--trackingVerbose $TrackingVerbosity"

  OtherOptionsGiven=true
else
  TrackingVerbosityString=""
fi


if [ ! "$TrackOpticalPhotons" == "" ] && ! $TrackOpticalPhotons ; then
  TrackOpticalPhotonsString="--noOpticalPhotonTracking"

  OtherOptionsGiven=true
else
  TrackOpticalPhotonsString=""
fi




# always use an init-file (in order to save the used program options)
if $InitFileExists && $OtherOptionsGiven ; then
  # abort the script
  echo "An existing init-file as well as other options are specified. The script is aborted"
  exit
# elif $InitFileExists && ! $OtherOptionsGiven ; then
  # do nothing, everything is OK
elif ! $InitFileExists && $OtherOptionsGiven ; then
  # write the options into the init-file
  if [ ! "$MacroFileString" == "" ] ; then
    echo $MacroFileString >> $InitFile
  fi
  if [ ! "$RunBatchModeString" == "" ] ; then
    echo $RunBatchModeString >> $InitFile
  fi
  if [ ! "$UsePhotonListString" == "" ] ; then
    echo $UsePhotonListString >> $InitFile
  fi
  if [ ! "$ParticleSourceMacroPathString" == "" ] ; then
    echo $ParticleSourceMacroPathString >> $InitFile
  fi
  if [ ! "$TileDimensionsString" == "" ] ; then
    echo $TileDimensionsString >> $InitFile
  fi
  if [ ! "$SearchOverlapsString" == "" ] ; then
    echo $SearchOverlapsString >> $InitFile
  fi
  if [ ! "$ScintillatorMaterialPropertiesString" == "" ] ; then
    echo $ScintillatorMaterialPropertiesString >> $InitFile
  fi
  if [ ! "$TwinTileReflectorMaterialPropertiesString" == "" ] ; then
    echo $TwinTileReflectorMaterialPropertiesString >> $InitFile
  fi
  if [ ! "$WrappingMaterialPropertiesString" == "" ] ; then
    echo $WrappingMaterialPropertiesString >> $InitFile
  fi
  if [ ! "$LGFibreMaterialPropertiesString" == "" ] ; then
    echo $LGFibreMaterialPropertiesString >> $InitFile
  fi
  if [ ! "$WLSFibreMaterialPropertiesString" == "" ] ; then
    echo $WLSFibreMaterialPropertiesString >> $InitFile
  fi
  if [ ! "$ScintiFibreMaterialPropertiesString" == "" ] ; then
    echo $ScintiFibreMaterialPropertiesString >> $InitFile
  fi
  if [ ! "$OpticalCementMaterialPropertiesString" == "" ] ; then
    echo $OpticalCementMaterialPropertiesString >> $InitFile
  fi
  if [ ! "$OutputDirectoryString" == "" ] ; then
    echo $OutputDirectoryString >> $InitFile
  fi
  if [ ! "$SimulationOutputFileNameString" == "" ] ; then
    echo $SimulationOutputFileNameString >> $InitFile
  fi
  if [ ! "$AddStringToOutputFilesString" == "" ] ; then
    echo $AddStringToOutputFilesString >> $InitFile
  fi
  if [ ! "$ControlVerbosityString" == "" ] ; then
    echo $ControlVerbosityString >> $InitFile
  fi
  if [ ! "$RunVerbosityString" == "" ] ; then
    echo $RunVerbosityString >> $InitFile
  fi
  if [ ! "$EventVerbosityString" == "" ] ; then
    echo $EventVerbosityString >> $InitFile
  fi
  if [ ! "$TrackingVerbosityString" == "" ] ; then
    echo $TrackingVerbosityString >> $InitFile
  fi
  if [ ! "$TrackOpticalPhotonsString" == "" ] ; then
    echo $TrackOpticalPhotonsString >> $InitFile
  fi

  echo >> $InitFile   # this is necessary to correctly read in the file later... Why? Don't ask me......
elif ! $InitFileExists && ! $OtherOptionsGiven ; then
  # abort the script
  echo "Neither an existing init-file nor other options are specified. The script is aborted"
  exit
fi




# finally, run the program
$BUILDDIR/RunSimulation $InitFileString
