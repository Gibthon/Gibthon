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
	var self = this;
	var url = '/fragment/get/' + id + '?value=meta';
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
	},
	_init: function(){
		console.log('ui.designer._init();');
		var self = this;
		this.canvas = this.element[0];
		this.ctx = this.canvas.getContext('2d');
		
		this.fragments = Array();
		
		/*Variables*/
		this.width = $(this.canvas).width();
		this.height = $(this.canvas).height();
		this.cx = this.width / 2;
		this.cy = this.height / 2;
		
		this.p_radius = Math.min(this.width, this.height) * 0.4;
		this.p_thickness = 20;
		this.p_margin = 5;
		
		this.length = 0;

		
		$(this.canvas).attr('height', this.height).attr('width', this.width);
		
		this._redraw();
	},
	addFragment: function(id){
		if(id == undefined) return;
		id = parseInt(id);
		if( isNaN(id) ) return;
		
		console.log('adding:' + id);
		var f = new Fragment(id);
		this.fragments.push(f);
		this.length = this.length + f.length;
		this._redraw();
	},
	_redraw: function(){
		this.ctx.clearRect(0,0,this.width, this.height);
		this._draw_title();	
		this._draw_plasmid();
	},
	_draw_title: function(){
		this.ctx.fillStyle = this.options.titleColour;
		this.ctx.font = this.options.titleFont;
		this.ctx.textAlign = 'center';
		
		var x = this.cx;
		var y = this.cy - 22;
		
		this.ctx.fillText(this.options.name, x, y);
		
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
			pos = pos + f.length * rad_per_bp;
		}
	},
	_draw_fragment: function(start, end, colour){
		console.log('_draw_fragment(' + start + ', ' + end + ')');
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
		console.log('angular_length: ' + angular_length);
		
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
