################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../core/borrar.cpp \
../core/flatfield.cpp \
../core/utility.cpp 

OBJS += \
./core/borrar.o \
./core/flatfield.o \
./core/utility.o 

CPP_DEPS += \
./core/borrar.d \
./core/flatfield.d \
./core/utility.d 


# Each subdirectory must supply rules for building sources it contributes
core/%.o: ../core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/cfitsio_lib_mac/include -I/Users/mrv/Documents/opencv-3.0.0/include/opencv -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


