# Porting alShaders to Arnold Render API v6 #
nominator for things to look at again:
> // ToDoJed: Fix for Porting->v6, 

state of the art and current progress in porting.

|               |       |                 |              |       |         |
|---------------|-------|-----------------|--------------|-------|---------| 
|shader         | ready | quality of port | improvements | C++11 | -Werror |
|---------------|-------|-----------------|--------------|-------|---------| 
|alSurface      |       |                 |              |       |         |             
|alRemap        |       |                 |              |       |         |           
|alCombine      |       |                 |              |       |         |             
|alNoise        |       |                 |              |       |         |           
|alLayer        |       |                 |              |       |         |           
|alPattern      |       |                 |              |       |         |             
|alCurvature    |       |                 |              |       |         |               
|alBlackbody    |       |                 |              |       |         |               
|alInputVector  |       |                 |              |       |         |                 
|alHair         |       |                 |              |       |         |          
|alColorSpace   |       |                 |              |       |         |                
|alTriplanar    |       |                 |              |       |         |               
|alSwitch       |       |                 |              |       |         |            
|alJitterColor  |       |                 |              |       |         |                 
|alCache        |       |                 |              |       |         |           
|alCel          |       |                 |              |       |         |         
|alFlake        |       |                 |              |       |         |           
|cryptomatte    |       |                 |              |       |         |               
|               |       |                 |              |       |         |    
|common         |       |                 |              |       |         |          
-----------------------------------------------------------------------------


## ToDo List ##

- make the dist folder result more descriptive  

       cmake-build-debug-visual-c-2015-x64/dist/1.0.0rc20/ai/bin

    better is with the arnold API version in the designation/folder:
    
       cmake-build-debug-visual-c-2015-x64/dist/alShaders-1.0.0rc20_Arnold-6.0.1.0/ai/bin
    
    This helps when using something like 
        
        ARNOLD_PLUGIN_PATH = "${ARNOLD_PLUGIN_PATH};D:/Users/Jedzia.pubsiX/htoa/alShaders-1.0.0rc20_Arnold-6.0.1.0/ai/bin"
    in Houdini.   
    


## Notes on Actions done
 
- what is #ifdef AI_GPU_COMPILER ? this disables AtString in some definitions?

- uses const AtString

        AiNodeGetPtr(options, "background"
        AiNodeIs(controller, AtString("uBasic_controller")))
        AiNodeGetBool(controller, AtString("strip_material_namespaces"));
        AiNodeGetBool(controller, AtString("override_cryptomatte"))) {
        AiNodeGetInt(controller, AtString("cryptomatte_depth"));
        AiNodeGetInt(controller, AtString("pointcloud_instance_verbosity"));
        AiAOVEnabled(


    alle 	const char- userName = AiShaderEvalParamStr(p_userName);
    zu 		AtString userName = AiShaderEvalParamStr(
    
    alle 	AiM4PointByMatrixMult(&result, mtx, &vector);
    zu		result = AiM4PointByMatrixMult(mtx, vector);	
    
    alle 	AiM4VectorByMatrixMult(&result, mtx, &vector);
    zu 		result = AiM4VectorByMatrixMult(mtx, vector);
    
    alle 	result = sg->area;
    zu 		result = AiShaderGlobalsArea(sg);
    
    all     std::string / .c_str() 
    to      AtString replacements 

- sg->Rr becomes sg->bounces
- sg->Rr_diff becomes sg->bounces_diffuse

- AiShaderEvalParamPnt becomes AiShaderEvalParamVec

- Replace	AiColorLerp() and AiRGBALerp() with AiLerp
- Replace	AiColorClamp()	with	the	templated	AiRGBClamp():	
- Replace	AiV3Exists()	with	AiV3IsFinite():		

- Replace	.c_str()	with	AtString replacements

## ToDo's about porting steps done

- Check everywhere `#if AI_VERSION_MAJOR_NUM > 0 || AI_VERSION_ARCH_NUM > 4` is if that 
  makes sense.
- ../common/alUtil.cpp was added to the alHair target because of AiEvaluateLightSample(...) proxy.
  Remove it, after cleanup.

## Hints about the theory ##

- PDF or Probability Density Function, see  3Blue1Brown, Why “probability of 0” does not 
  mean “impossible” | Probabilities of probabilities, part 2 - 
  https://www.youtube.com/watch?v=ZA4JkHKZM50 

## Hints for compilation
#### Linux
    export ARNOLD_PATH=/path/to/arnold
    c++ simple_shader.cpp -o simple_shader.so -Wall -O2 -shared -fPIC -I$ARNOLD_PATH/include -L$ARNOLD_PATH/bin -lai
#### macOS
    export ARNOLD_PATH=/path/to/arnold
    c++ simple_shader.cpp -o simple_shader.dylib -Wall -O2 -shared -fPIC -I$ARNOLD_PATH/include -L$ARNOLD_PATH/bin -lai
#### Windows Visual Studio command prompt
    set ARNOLD_PATH=c:/path/to/arnold
    cl /LD simple_shader.cpp /I %ARNOLD_PATH%/include %ARNOLD_PATH%/lib/ai.lib /link /out:simple_shader.dll