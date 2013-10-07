import uigen

ui.shader({
	'name':'alRemapFloat',
	'description':'Controls to adjust a float value',
	'output':'rgb',
	'maya_name':'alRemapFloat',
	'maya_classification':'texture',
	'maya_id':'0x0011640C',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_RemapColor',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('input', 'float', 'Input')

uigen.remapControls(ui)