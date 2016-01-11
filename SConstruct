import os
import sys

###############################################
# Options
###############################################

AddOption('--verbose', action='store_true', dest='verbose', default=False, help='Enable print statements')
AddOption('--cuda', dest='cuda', help='Compile gpu components')
AddOption('--openmp', action='store_true', dest='openmp', help='Compile openmp components')

###############################################
# Core
###############################################

VariantDir('build/core', 'core')

envCore = Environment()
envCore.Append(CCFLAGS=' -Ofast')
envCore.Append(CCFLAGS=' -std=c99')
envCore.Replace(LIBS=[])

envCore.Append(CCFLAGS=' -DSILENCE='+str(0 if GetOption('verbose') else 1))

# Library compile defines/paths

if GetOption('openmp'):
	envCore.Append(CCFLAGS=' -DOPENMP')
	envCore.Append(CCFLAGS=' -fopenmp')
	envCore.Append(LINKFLAGS=' -fopenmp')

# Platform-dependent build
# (may require library-specific changes)

objects = []
if sys.platform == 'darwin':
	envCore.Append(CCFLAGS=' -DMAC_OS')
	envCore.Replace(CC='clang')
	envCore.Append(CCFLAGS=' -Wno-incompatible-pointer-types')
	envCore.Append(LINKFLAGS=' -framework Accelerate')

	objects.append(envCore.SharedObject('build/core/paradmm-host', Glob('build/core/*.c')))

else:
	envCore.Replace(CC='gfortran')
	envCore.Append(CCFLAGS=' -Wno-implicit')
	envCore.Append(CCFLAGS=' -Wno-unused-result')

	objects.append(envCore.SharedObject('build/core/paradmm-host', Glob('build/core/*.c')))
	
envCore.SharedLibrary('out/lib/paradmm', objects)

###############################################
# Cuda
###############################################

envNVCC = None
o = None
if GetOption('cuda'):
	envNVCC = Environment()

	envNVCC['CUDA_TOOLKIT_PATH'] = GetOption('cuda')
	envNVCC['CUDA_SDK_PATH'] = GetOption('cuda')
	envNVCC.Tool('cuda','.')
	
	if GetOption('openmp'):
		envNVCC.Append(NVCCFLAGS=' -DOPENMP')
		envNVCC.Append(NVCCFLAGS=' -Xcompiler -fopenmp')
		envNVCC.Append(LINKFLAGS=' -Xcompiler -fopenmp')
		
	envNVCC.Append(NVCCFLAGS=' -Xcompiler -Ofast')
	envNVCC.Append(NVCCFLAGS=' -Xcompiler -fPIC')
	envNVCC.Append(NVCCFLAGS=' --use_fast_math -m64 --shared -dc')
	envNVCC.Append(NVCCFLAGS=' -DSILENCE='+str(0 if GetOption('verbose') else 1))
	envNVCC.Append(NVCCFLAGS=' -Ibuild/core')
	envNVCC.Append(NVCCFLAGS=' -DNEEDS_EXTERN_C')
	envNVCC.Replace(CC='nvcc')
	
	envNVCC['STATIC_AND_SHARED_OBJECTS_ARE_THE_SAME']=1
	o = envNVCC.Object('out/lib/paradmm-gpu', Glob('build/core/*.cu'))
	
	# envNVCC.SharedLibrary('out/lib/paradmm-gpu', o)

###############################################
# Includes
###############################################

envCore.Install('out/include', Glob('build/core/*.h'))

###############################################
# Packing Demo
###############################################

VariantDir('build/packing', 'packing')

envPacking = envCore.Clone()
envPacking.Append(CPPPATH=['#out/include'])
envPacking.Append(LIBPATH=['#out/lib'])
envPacking.Append(LIBS='paradmm')

envPacking.Program('#out/packing', Glob('build/packing/*.c'))

##

if envNVCC is not None:
	envPackingNVCC = envNVCC.Clone()
	
	envPackingNVCC.Append(LIBPATH=['#out/lib'])
	# envPackingNVCC.Append(LIBS=['paradmm', 'paradmm-gpu'])
	envPackingNVCC.Append(LIBS=['paradmm'])
	
	o2 = envPackingNVCC.Object(Glob('build/packing/*.cu'))
	# envPackingNVCC.Program('#out/packing-gpu', [o2])
	envPackingNVCC.Program('#out/packing-gpu', [o, o2])
