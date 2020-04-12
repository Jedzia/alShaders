# Porting alShaders to Arnold Render API v6 #
nominator for things to look at again:
> ToDoJed: Fix for Porting->v6, 

# Notes on Actions done
 
- what is #ifdef AI_GPU_COMPILER ? this disables AtString in some definitions?

- uses const AtString

        AiNodeGetPtr(options, "background"
        AiNodeIs(controller, AtString("uBasic_controller")))
        AiNodeGetBool(controller, AtString("strip_material_namespaces"));
        AiNodeGetBool(controller, AtString("override_cryptomatte"))) {
        AiNodeGetInt(controller, AtString("cryptomatte_depth"));
        AiNodeGetInt(controller, AtString("pointcloud_instance_verbosity"));
        AiAOVEnabled(


    alle 	const char* userName = AiShaderEvalParamStr(p_userName);
    zu 		AtString userName = AiShaderEvalParamStr(
    
    alle 	AiM4PointByMatrixMult(&result, mtx, &vector);
    zu		result = AiM4PointByMatrixMult(mtx, vector);	
    
    alle 	AiM4VectorByMatrixMult(&result, mtx, &vector);
    zu 		result = AiM4VectorByMatrixMult(mtx, vector);
    
    alle 	result = sg->area;
    zu 		result = AiShaderGlobalsArea(sg);

* sg->Rr becomes sg->bounces

* AiShaderEvalParamPnt becomes AiShaderEvalParamVec

* Replace	AiColorLerp() and AiRGBALerp() with AiLerp
* Replace	AiColorClamp()	with	the	templated	AiRGBClamp():	
* Replace	AiV3Exists()	with	AiV3IsFinite():		


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