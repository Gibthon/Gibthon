(function( $, undefined ) {

var colnum = 0;

var hsv_2_rgb = function(h,s,v)
{
	var r, g, b;

	if(s==0)
	{
		r=g=b=Math.round(v*255);
	}
	else
	{
		// h must be < 1
		var var_h = h * 6;
		if (var_h==6) var_h = 0;
		//Or ... var_i = floor( var_h )
		var var_i = Math.floor( var_h );
		var var_1 = v*(1-s);
		var var_2 = v*(1-s*(var_h-var_i));
		var var_3 = v*(1-s*(1-(var_h-var_i)));
		if(var_i==0)
		{ 
			var_r = v; 
			var_g = var_3; 
			var_b = var_1;
		}
		else if(var_i==1)
		{ 
			var_r = var_2;
			var_g = v;
			var_b = var_1;
		}
		else if(var_i==2)
		{
			var_r = var_1;
			var_g = v;
			var_b = var_3
		}
		else if(var_i==3)
		{
			var_r = var_1;
			var_g = var_2;
			var_b = v;
		}
		else if (var_i==4)
		{
			var_r = var_3;
			var_g = var_1;
			var_b = v;
		}
		else
		{ 
			var_r = v;
			var_g = var_1;
			var_b = var_2
		}
		//rgb results = 0 ÷ 255  
		r=Math.round(var_r * 255);
		g=Math.round(var_g * 255);
		b=Math.round(var_b * 255);
	}
	return {r:r, g:g, b:b}
}

var get_col = function(){
	var h = (colnum * (5 / 7)) % 1;
	var s = 0.7;
	var v = 1.0;
	
	var col = hsv_2_rgb(h,s,v);
	
	colnum = colnum + 1;
	return 'rgb('+col.r+','+col.g+','+col.b+')';
}

var PI2 = 2 * Math.PI;


$.widget("ui.designer", {
	options: {
		name: '',
		id: 0,
		titleFont: '20px bold Arial',
		titleColour: 'rgb(20,20,20)',
		lengthFont: '16px Arial',
		lengthColour: 'rgb(100,100,100)',
	},
	_init: function(){
	},
	_create: function(){		
		var self = this;
		this.canvas = this.element[0];
	},
// ------------------------------------------------------------------------------------ PUBLIC API
	addFragment: function(id, tell_server, redraw){	
	},
	changeName: function(new_name){
	},
// ---------------------------------------------------------------------------------- PRIVATE API

});

})(jQuery);
