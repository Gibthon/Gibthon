(function( $, undefined ) {

var STATE_NORMAL = 0;
var STATE_REORDER = 0;
var STATE_ROTATE = 0;

var cols = [
	'rgb(255,0,0)',
	'rgb(0,255,0)',
	'rgb(0,0,255)',
	'rgb(255,255,0)',
	'rgb(255,0,255)',
	'rgb(0,255,255)',
	];

var get_col = function(num){
	return cols[ num % cols.length ];
}

function Fragment(id)
{
	this.id = id;
	this.center = 0; //the angle of the center of the fragment
	var self = this;
	var url = '/fragment/get/' + id + '/?value=meta';
	$.ajax({
	  url: url,
	  dataType: 'json',
	  async: false,
	  success: function(data)
	  {
		  self.name = data[1].name;
		  self.desc = data[1].desc;
		  self.length = data[1].length;
	  }
	});
}

$.widget("ui.designer", {
	options: {
		name: '',
		id: 0,
		titleFont: '20px bold Arial',
		lengthFont: '16px Arial',
		titleColour: 'rgb(20,20,20)',
		lengthColour: 'rgb(100,100,100)',
		labelAreaWidth: 0.25, 
		labelFont: '16px italic Arial',
		lineColour: 'rgb(100,100,100)',
		stubLength: 20,
	},
	_init: function(){
		console.log('ui.designer._init();');
		this._redraw();
	},
	_create: function(){
		console.log('ui.designer._create();');
		var self = this;
		this.canvas = this.element[0];
		this.ctx = this.canvas.getContext('2d');
		
		this.fragments = Array();
		
		/*Create Layout Variables*/
		this.width = 0;
		this.height = 0;
		this.cx = 0;
		this.cy = 0;
		
		this.p_radius = 0;
		this.p_thickness = 0;
		
		this.la_width = 0;//label area width
		this.ll_up = 0;//max height for label-line to be upwards
		this.ll_down = 0;//min height for label-line to be downwards
		
		//Calculate what they should be.
		this._update_layout();
		
		
		this.length = 0;
		this.name = this.options.name;

		$.getJSON('fragments/', [], function(data) {
			if(data[0] != 0)
			{
				console.log('could not get initial fragments: ' + data[1]);
				return;
			}
			console.log('adding initial fragments');
			for(var i in data[1])
			{
				self.addFragment(data[1][i], false, false);
			}
			self._redraw();
		});


		this._redraw();
	},
	addFragment: function(id, tell_server, redraw){
		if(id == undefined) return;
		if(tell_server == undefined) tell_server = true;
		if(redraw == undefined) redraw = true;
		id = parseInt(id);
		if( isNaN(id) ) return;
		
		var f = new Fragment(id);
		this.fragments.push(f);
		this.length = this.length + f.length;
		
		if(tell_server)
		{
			$.getJSON('addFragment/' + id + '/', [], function(data) {
				if(data[0] != 0)
					console.log('Error while adding fragment: ' + data[1]);
			});
		}
		if(redraw)
			this._redraw();
	},
	changeName: function(new_name){
		this.name = new_name;
		this._redraw();
	},
	_update_layout: function() {
		this.width = $(this.canvas).width();
		this.height = $(this.canvas).height();
		$(this.canvas).attr('height', this.height).attr('width', this.width);
		
		this.cx = this.width / 2;
		this.cy = this.height / 2;
		
		this.p_radius = Math.min(this.width * ( 1 - 2 * this.options.labelAreaWidth), this.height) * 0.4;
		this.p_thickness = this.p_radius * 0.1;
		
		this.la_width = this.width * this.options.labelAreaWidth;
		this.ll_up = this.cy - this.p_radius * 0.7071;
		this.ll_down = this.cy + this.p_radius * 0.7071;
	},
	_redraw: function(){
		this.ctx.clearRect(0,0,this.width, this.height);
		this._draw_title();	
		this._draw_plasmid();
		this._draw_labels();
	},
	_draw_title: function(){
		this.ctx.fillStyle = this.options.titleColour;
		this.ctx.font = this.options.titleFont;
		this.ctx.textAlign = 'center';
		
		var x = this.cx;
		var y = this.cy - 22;
		
		this.ctx.fillText(this.name, x, y);
		
		this.ctx.fillStyle = this.options.lengthColour;
		this.ctx.font = this.options.lengthFont;
		
		var len_string = this.length + ' bp';
		var x = this.cx;
		var y = this.cy + 2;
		
		this.ctx.fillText(len_string, x, y);
	},
	_draw_plasmid: function(){
		if(this.fragments.length == 0)
		{
			this._draw_blank(0, 2 * Math.PI);
			return;
		}
		var pos = 0;
		var rad_per_bp = 2 * Math.PI / (this.length);
		
		for(var i  in this.fragments)
		{
			var f = this.fragments[i];
			this._draw_fragment(pos, pos + f.length * rad_per_bp, get_col(i));
			f.center = pos + (f.length * rad_per_bp / 2);
			pos = pos + f.length * rad_per_bp;
		}
	},
	_draw_labels: function(){
		var left = this.la_width / 2;
		var right = this.width - this.la_width / 2;
		
		this.ctx.fillStyle = this.options.titleColour;
		this.ctx.font = this.options.labelFont;
		this.ctx.textAlign = 'center';
		
		for(var i  in this.fragments)
		{
			var f = this.fragments[i];
			var sy = this.cy + this.p_radius * Math.sin(f.center);
			var sx = this.cx + this.p_radius * Math.cos(f.center);
			var ex = 0; var ey = 0;
			if(sx < this.cx) //label on left
			{
				//horizontal position
				ex = left + (2 + 0.5 * this.ctx.measureText(f.name).width);
				var text_x = left;
				//vertical position
				ey = sy;
			}
			else //label on right
			{
				//horizontal position
				ex = right - (2 + 0.5 * this.ctx.measureText(f.name).width);
				var text_x = right;
				//vertical position
				ey = sy;
			}
			var text_y = ey;
			this._draw_handle(sx, sy);
			
			this.ctx.save();
			
			this.ctx.strokeStyle = this.options.lineColour;
			this.ctx.lineWidth = 1;
			
			this.ctx.beginPath();			
			//draw stub
			if(sy < this.ll_up) //stub is upwards
			{
				this.ctx.moveTo(sx, sy - 5);
				text_y = text_y - this.options.stubLength;
				sy = sy - this.options.stubLength;
			}
			else if(sy > this.ll_down) //stub is downwards
			{
				this.ctx.moveTo(sx, sy + 5);
				text_y = text_y + this.options.stubLength;
				sy = sy + this.options.stubLength;
			}
			else if(sx > this.cx) //stub is right
			{
				this.ctx.moveTo(sx + 5, sy);
				sx = sx + this.options.stubLength;
			}
			else if(sx < this.cx) //stub is left
			{
				this.ctx.moveTo(sx - 5, sy);
				sx = sx - this.options.stubLength;
			}
			
			this.ctx.lineTo(sx, sy);
			
			//draw to text
			
			this.ctx.lineTo(ex,sy);
			
			
			this.ctx.stroke();
			
			this.ctx.restore();
			
			this.ctx.fillText(f.name, text_x, text_y + 6);
		}
		
	},
	_draw_handle: function(x, y){
		this.ctx.save();
		this.ctx.fillStyle = this.options.lineColour;
		this.ctx.strokeStyle = this.options.lineColour;
		this.ctx.lineWidth = 1;
		
		this.ctx.beginPath();
		this.ctx.arc(x,y,3,0,2*Math.PI);
		this.ctx.fill();
		
		this.ctx.beginPath();
		this.ctx.arc(x,y,5,0,2*Math.PI);
		this.ctx.stroke();
		
		this.ctx.restore();
	},
	_draw_fragment: function(start, end, colour){
		this.ctx.save();
		this.ctx.strokeStyle = colour;
		this.ctx.lineWidth = this.p_thickness;
		var space = 1.5 / this.p_radius;
		
		if(start == end) end = end + 2 * Math.PI;
		
		this.ctx.beginPath();
		
		this.ctx.arc(this.cx, this.cy, this.p_radius, start + space, end - space);
		
		this.ctx.stroke();
		
		this.ctx.restore();
	},
	_draw_blank: function(start, end){
		this.ctx.save();
		this.ctx.strokeStyle = 'rgb(100,100,100)';
		if(start == end) end = end + 2 * Math.PI;
		this._dotted_arc(this.cx, this.cy, this.p_radius - this.p_thickness / 2, start, end);
		this._dotted_arc(this.cx, this.cy, this.p_radius + this.p_thickness / 2, start, end);
		this.ctx.restore();
	},
	_dotted_arc: function(x, y, radius, start, end, angular_length){
		if((angular_length == undefined) || (angular_length <= 0))
			angular_length = Math.PI / 50;
		var p = start;
		
		while((p + angular_length) < end)
		{
			this.ctx.beginPath();
			this.ctx.arc(x,y,radius, p, p + angular_length);
			this.ctx.stroke();
			p = p + 2 * angular_length;
		}
		if((end - p) < angular_length)
			this.ctx.arc(x,y,radius, p, end);
	}
	
});

})(jQuery);
