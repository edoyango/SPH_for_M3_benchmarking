@echo off
	
set src_dir=.\src\
set obj_dir=.\obj\

if not exist %obj_dir% mkdir %obj_dir%

REM checking if debug is passed to script if so, compile with /check option
set compoption=/traceback /O3 /Qxhost /Qipo /Qopenmp
:loop
	if "%1"=="debug" ( 
		set compoption=/traceback /O3 /Qxhost /Qipo /check:all /Qopenmp )
	if "%1"=="clean" (
		goto clean )
	if "%1"=="defaultoption" (
		set compoption=/traceback )
shift
if not "%~1"=="" goto loop

mpiifort ^
/object:%obj_dir% /module:%obj_dir% ^
%compoption% ^
%src_dir%constants.f90 ^
%src_dir%param.f90 ^
%src_dir%datatypes.f90 ^
%src_dir%globvar.f90 ^
%src_dir%globvar_para.f90 ^
%src_dir%param_para.f90 ^
%src_dir%summary_m.f90 ^
%src_dir%error_msg_m.f90 ^
%src_dir%material_rates_m.f90 ^
%src_dir%output_m.f90 ^
%src_dir%input_m.f90 ^
%src_dir%kernel_m.f90 ^
%src_dir%flink_lisk_m.f90 ^
%src_dir%ORB_sr_m.f90 ^
%src_dir%ORB_m.f90 ^
%src_dir%single_step_m.f90 ^
%src_dir%time_integration_m.f90 ^
%src_dir%main.f90 ^
/exe:.\sph.exe

goto:eof

:clean
del %obj_dir%\*.obj
del %obj_dir%\*.mod
del %obj_dir%\*.o