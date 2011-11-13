(function( $, undefined ) {

var STATE_NORMAL = 0;
var STATE_REORDER = 1;
var STATE_ROTATE = 2;

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

/*
 * Determine whether angle is between lower and higher.
 * if below lower, return -1,
 * if within return 0
 * if above higher return 1
 * */
var is_in = function(lower, higher, angle)
{
	//relative to lower
	var h = (higher - lower);
	var a = (angle - lower);
	//in range 0 - 2*PI
	while(h < 0) h = h + PI2;
	while(h > PI2) h = h - PI2;
	while(a < 0) a = a + PI2;
	while(a > PI2) a = a - PI2;
	
	if(a < h) return 0;
	
	if( (a - h) > (PI2 - a) )
		return -1; //closer to lower
	return 1; //closer to higher
}

/*
 * Contain all the information needed for fragment display
 * */
function Fragment(data)
{	
	var self = this;
	this.id = data.fid;
	this.cfid = data.cfid;
	this.name = data.name;
	this.desc = data.desc;
	this.length = data.length;
		
	/* Layout / Postitioning based */
	this.start = 0; //the start angle of the fragment
	this.end = 0;
	this.center = function() { return (self.start + self.end) / 2.0; };
	this.color = get_col();	
	
	this.highlight = false;
	
	this.reordering = false;
	this.handle_angle = 0; //the offset from the center of the fragment where the mouse originally clicked
	this.reorder_pos = 0;
	return this;
}

/*
 * Fetch info on a fragment from the server
 * */
function get_fragment_from_id(id)
{
	var url = '/fragment/get/' + id + '/?value=meta';
	var f;
	$.ajax({
		url: url,
		dataType: 'json',
		async: false,
		success: function(data)
		{
			if(data[0] != 0)
			{
				throw 'Error while getting fragment data: ' + data[1];
			}
			f = new Fragment(data[1]);
			for(var z in f)
			return f;
		}
	});
	return f;
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
		clickTime: 250,
	},
	_init: function(){
		//this._redraw();
	},
	_create: function(){
		this.state = STATE_NORMAL;
		
		var self = this;
		this.canvas = this.element[0];
		this.$canvas = $(this.canvas);
		this.ctx = this.canvas.getContext('2d');
		
		this.fragments = Array();
		
		/*Create Layout Variables*/
		this.width = 0;
		this.height = 0;
		this.cx = 0;
		this.cy = 0;
		
		this.p_radius = 0;
		this.p_thickness = 0;
		this.outer_radius = 0;
		
		this.la_width = 0;//label area width
		this.ll_up = 0;//max height for label-line to be upwards
		this.ll_down = 0;//min height for label-line to be downwards
		
		//Calculate what they should be.
		this._update_layout();
		
		
		this.length = 0;
		this.name = this.options.name;
		this.rad_per_bp = 0;
		this.next = 0; this.prev = 0;

		$.getJSON('fragments/', [], function(data) {
			if(data[0] != 0)
			{
				console.log('could not get initial fragments: ' + data[1]);
				return;
			}

			for(var i in data[1])
			{
				self._add_fragment(new Fragment(data[1][i]));
			}
			self._layout_fragments();
			self._redraw();
		});

		this._update = true; //should the whole thing be redrawn?
		
		/* Mouse things*/
		this.$canvas.mousemove( function(event) { self._mouse_move(event); });
		this.$canvas.mousedown( function(event) { self._mouse_down(event); });
		this.$canvas.mouseup( function(event) { self._mouse_up(event); });
		this.$canvas.mouseleave( function(event) { self._mouse_up(event); });
		
		this._mdown = false;
		this._mdown_time = 0;
		this._mdown_pos = {r:0, theta: 0,};
		this._mouse_over = -1; //which fragment is the mouse over? -1 = none
		
		/* add the clipping dialog*/
		this.$clipping = this.$canvas.siblings('#fragment_clipping');
		this.$clipping.dialog({ autoOpen: false, title: 'Clipping', modal: true, resizable: false});
		
		/*Add the info box*/
		this.$info = this.$canvas.siblings('.fragment-info');
		this.$info.find('#fragment_clip').button({
			icons: {primary: 'ui-icon-scissors'},
		}).click( function() {
			var fid = parseInt(self.$info.attr('fragment'));
			var cfid = self.fragments[fid].cfid;
			console.log('load("clipping/' + cfid + '/", function(responseText, textStatus, XMLHttpRequest) {...});');
			self.$clipping.load('clipping/' + cfid + '/', function() {
				console.log('hello!');
				self.$clipping.dialog({'width': Math.max($(window).width() * 0.5, 500)});
				self.$clipping.dialog('open');
			});
			self.$info.hide();
		});
		this.$info.find('#fragment_remove').button({
			icons: {primary: 'ui-icon-trash'},
		}).click( function() {
			var f = parseInt(self.$info.attr('fragment'));
			self._remove_fragment(f);
			self.$info.hide();
			self._redraw();
		});
		
		
	},
// ------------------------------------------------------------------------------------ PUBLIC API
	addFragment: function(id, tell_server, redraw){
		var self = this;
		if(id == undefined) return;
		if(tell_server == undefined) tell_server = true;
		if(redraw == undefined) redraw = true;
		id = parseInt(id);
		if( isNaN(id) ) return;
		
		var f = get_fragment_from_id(id);
		
		var fid = this._add_fragment(f);

		if(tell_server)
		{
			$.getJSON('addFragment/' + id + '/', [], function(data) {
				if(data[0] != 0)
					throw ('Error while adding fragment: ' + data[1]);
				else
				{
					self.fragments[fid].cfid = data[1].cfid;
				}
			});
		}
		
		if(redraw)
			this._redraw();
	},
	changeName: function(new_name){
		this.name = new_name;
		this._update = true;
		this._redraw();
	},
// ---------------------------------------------------------------------------------- PRIVATE API
/*
 *  ----------------------------- Fragment manipulation -----------------------------
 * */
	_add_fragment: function(fragment)
	{
		if(fragment != undefined)
		{
			this.fragments.push(fragment);
			this.length = this.length + fragment.length;
			this.rad_per_bp = 2 * Math.PI / (this.length);
			this._update = true;
			this._layout_fragments();
			return (this.fragments.length -1);
		}
		return -1;
	},
	_remove_fragment: function(fid){
		if(fid < 0) return;
		this.length = this.length - this.fragments[fid].length;
		this.rad_per_bp = PI2 / (this.length);
		
		$.getJSON('removeFragment/' + this.fragments[fid].cfid + '/', function(data) {
			if(data[0] != 0)
				alert('Error while removing fragment: ' + data[1]);
		});
		this.fragments.splice(fid, 1);
		this._update = true;
		this._layout_fragments();
	},
	_swap_fragments: function(a, b) //assume that a is the fragment before b
	{
		var fa = this.fragments[a];
		var fb = this.fragments[b];
		
		fb.start = fa.start;
		fb.end = fb.start + fb.length * this.rad_per_bp;
		fa.start = fb.end;
		fa.end = fa.start + fa.length * this.rad_per_bp;
				
		this.fragments[a] = fb;
		this.fragments[b] = fa;
		
		this._redraw();
	},
	_update_order: function(){ //let the server know that the order changed
		var order = new Array();
		for(i in this.fragments)
		{
			order.push(this.fragments[i].cfid);
		}
		$.getJSON('saveOrder/', {'order[]': order,}, function(data) {
			$('#status').text(data[1].modified);
		});
	},
/*
 * ----------------------------- Drawing -----------------------------
 * */
	_redraw: function(){
		if(this._update)
		{
			this.ctx.clearRect(0,0,this.width, this.height);
			this._draw_title();	
			this._draw_plasmid();
			
			this._draw_labels();
				
			this._update = false;
		}
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
				
		for(var i  in this.fragments)
		{
			this._draw_fragment(i);
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
			var sy = this.cy + this.p_radius * Math.sin(f.center());
			var sx = this.cx + this.p_radius * Math.cos(f.center());
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
	_draw_fragment: function(f_id){	
		var f = this.fragments[f_id];
		this.ctx.save();
		this.ctx.strokeStyle = f.color;
		this.ctx.lineWidth = this.p_thickness;
		if(f.highlight)
		{
			this.ctx.globalAlpha = 0.5;
		}
				
		var r = this.p_radius;
		var s = f.start;
		var e = f.end;

		if(f.reordering)
		{
			this._draw_blank(s, e);
			r = this.outer_radius;
			s = f.reorder_pos - (f.length / 2 * this.rad_per_bp);
			e = f.reorder_pos + (f.length / 2 * this.rad_per_bp);
		}
		
		this.ctx.beginPath();
		this.ctx.arc(this.cx, this.cy, r, s, e);
		this.ctx.stroke();
		this.ctx.restore();
	},
	_draw_blank: function(start, end){
		this.ctx.save();
		this.ctx.strokeStyle = 'rgb(100,100,100)';
		this.ctx.lineWidth = 1;
		this.ctx.globalAlpha = 1.0;
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
	},
/*
 * ----------------------------- Layout -----------------------------
 * */
	_update_layout: function() {
		this.width = $(this.canvas).width();
		this.height = $(this.canvas).height();
		$(this.canvas).attr('height', this.height).attr('width', this.width);
		
		this.cx = this.width / 2;
		this.cy = this.height / 2;
		
		this.p_radius = Math.min(this.width * ( 1 - 2 * this.options.labelAreaWidth), this.height) * 0.4;
		this.p_thickness = this.p_radius * 0.1;
		this.outer_radius = this.p_radius + 1.5 * this.p_thickness;
		
		this.la_width = this.width * this.options.labelAreaWidth;
		this.ll_up = this.cy - this.p_radius * 0.7071;
		this.ll_down = this.cy + this.p_radius * 0.7071;
	},
	_layout_fragments: function() //layout the fragments in their order
	{
		if(this.fragments.length < 1) return;
		var pos = this.fragments[0].start;
		
		for(var f_id in this.fragments)
		{
			var f = this.fragments[f_id];
			f.start = pos * this.rad_per_bp;
			pos = pos + f.length;
			f.end = pos * this.rad_per_bp;
		}
	},
/*
 * ----------------------------- INTERACTION CALLBACKS -----------------------------
 * */
	_fragment_enter: function(f, pos) { //called when the mouse enters a fragment
		this.fragments[f].highlight = true;
		this._set_cursor('pointer');
		this._update = true;
	},
	_fragment_leave: function(f, pos) { //called when the mouse leaves a fragment
		this.fragments[f].highlight = false;
		this._set_cursor();
		this._update = true;
	},
	_fragment_click: function(f, pos) { //called when the user clicks on a fragment
		this.$info.find('#fragment_name').text(this.fragments[f].name);
		this.$info.find('#fragment_desc').text(this.fragments[f].desc);
		this.$info.attr('fragment', f);
		
		var ang = this.fragments[f].center();
		var x = this.cx + this.p_radius * Math.cos(ang) - 35;
		var y = this.cy + this.p_radius * Math.sin(ang) - (this.$info.height() + 16);
		
		this.$info.css({'left': x, 'top':y});
		this.$info.show();
	},
	_click: function(pos) { //called when a region which isn't a fragment is clicked
		this.$info.hide();
	},
	_fragment_init_drag: function(sel, pos) {
		this.$info.hide();
		this.state = STATE_REORDER;
		var f = this.fragments[sel];
		f.highlight = false;
		f.reordering = true;
		f.handle_angle = this._mdown_pos.theta - f.center();
		f.reorder_pos = f.center();
		this._set_cursor('move');
		this._updated = true;
	},
	_fragment_drag: function(pos) {
		for(i in this.fragments)
		{
			i = parseInt(i);
			var f = this.fragments[i];
			if(f.reordering)
			{
				f.reorder_pos = pos.theta - f.handle_angle;
				
				//check for reordering
				var n = i + 1;
				var p = i - 1;
				if(n >= this.fragments.length) n = 0;
				if(p < 0) p = this.fragments.length - 1;
				
				var next = this.fragments[n].center();
				var prev = this.fragments[p].center();
				
				switch( is_in(prev, next, pos.theta) )
				{
					case -1:
						this._swap_fragments(p, i);
						break;
					case 1:
						this._swap_fragments(i, n);
						break;
					default:
						break;
				}
				
				this._update = true;
			}
		}
	},
	_fragment_end_drag: function(f, pos) {
		this._update_order();
		this.state = STATE_NORMAL;
		this._set_cursor();
		for(i in this.fragments)
			this.fragments[i].reordering = false;
		this._updated = true;
	},
/* 
 * ----------------------------- MOUSE FUNCTIONS -----------------------------
 * */
	_get_xy: function(event) {
		var offset = this.$canvas.offset();
		return {x: event.pageX - offset.left, y: event.pageY - offset.top};
	},
	_get_rtheta: function(event) {
		var xy = this._get_xy(event);
		var x = xy.x - this.cx;
		var y = xy.y - this.cy;
		var t = Math.atan2(y, x);
		if( t < 0) t = 2 * Math.PI + t;
		return {r: Math.sqrt(x*x + y*y), theta: t};
	},
	_set_cursor: function(cursor) {
		if(cursor == undefined) cursor = 'default';
		document.body.style.cursor = cursor;
	},
	_get_selected_fragment: function(theta) { //given theta, find which fragment is moused over
		for(var i in this.fragments)
		{
			var f = this.fragments[i];
			if( is_in(f.start, f.end, theta) == 0)
			{
				return parseInt(i);
			}
		}
		return -1;
	},
	_is_in_plasmid: function(r) //is the mouse within the fragment?
	{
		return (r > (this.p_radius - this.p_thickness/2)) && (r < (this.p_radius + this.p_thickness/2));
	},
	_mouse_move: function(event){
		var pos = this._get_rtheta(event);
		switch(this.state)
		{
			case STATE_NORMAL:
				if(this._is_in_plasmid(pos.r))
				{
					if(!this._mdown)
					{
						//Mouse is not down and the mouse is in a fragment
						var f = this._get_selected_fragment(pos.theta);
						if( f < 0 ) throw ('_get_selected_fragment(' + pos.theta + ') returned "' + f + '"');
						if( this._mouse_over != f )
						{
							if(this._mouse_over >= 0) this._fragment_leave(this._mouse_over);
							this._mouse_over = f;
							if(f >= 0) this._fragment_enter(f);
						}
					}
					else
					{
						//if this is too long to be a click
						if((event.timeStamp - this._mdown_time) > this.options.clickTime)
						{
							//start reordering
							var selected = this._get_selected_fragment(this._mdown_pos.theta);
							if(selected >= 0)
							{
								this._fragment_init_drag(selected);
							}
							else
								throw ('Error initting reordering: selected = ' + selected);
						}
					}
				}
				else if(this._mouse_over >= 0)
				{
					this._fragment_leave(this._mouse_over);
					this._mouse_over = -1;
				}
				break;
			
			case STATE_REORDER:
				//handle dragging
				this._fragment_drag(pos);
				break;
		}
		this._redraw();
	},
	_mouse_down: function(event){
		this._mdown_pos = this._get_rtheta(event);
		this._mdown_time = event.timeStamp;
		this._mdown = true;
	},
	_mouse_up: function(event){
		var pos = this._get_rtheta(event);
		this._mdown = false;
		if((event.timeStamp - this._mdown_time) < this.options.clickTime)
		{
			if(this._is_in_plasmid(pos.r))
			{
				var i = this._get_selected_fragment(pos.theta);
				if( i < 0 ) throw 'this._get_selected_fragment(' + pos.theta + ') returned ' + i;
				this._fragment_click(i, pos);
				this._redraw();
			}
			else
			{
				this._click(pos);
				this._redraw();
			}
		}
		else if(this.state != STATE_NORMAL)
		{
			this._fragment_end_drag();
			this._redraw();
		}
		
	},
	
});

})(jQuery);
