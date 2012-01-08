var canvas;
var stage;


//Declare layout values
//constants
var LIB_HANDLE_WIDTH = 20;
var LIB_STOPPER_WIDTH = 6;
var LIB_ITEM_HEIGHT = 20;
var LIB_ITEM_SPACING = 4;

var FRAG_W = 25;
var FRAG_W_2 = FRAG_W / 2.0;
var FRAG_ARROW = 5.0;

//variable
var bounds = new Rectangle();
var LIB_X = 0; var LIB_Y = 0; var LIB_W = 0; var LIB_H; var LIB_SLIDE = 0;
var CW_R = 0; var CCW_R = 0; var DRAG_R = 30;

//colours
var UI_BG_FILL 				= Graphics.getRGB(108,165,208);
//var UI_BG_FILL_HL			= Graphics.getRGB();
//var UI_BG_FILL_SL			= Graphics.getRGB();

var UI_BG_STROKE 			= Graphics.getRGB(66,151,215);
//var UI_BG_STROKE_HL			= Graphics.getRGB();
//var UI_BG_STROKE_SL			= Graphics.getRGB();

var UI_FG_FILL 				= Graphics.getRGB(230,242,252);
var UI_FG_FILL_HL 			= Graphics.getRGB(208,229,245);
var UI_FG_FILL_SL		 	= Graphics.getRGB(246,248,249);

var UI_FG_STROKE 			= Graphics.getRGB(197,219,236);
var UI_FG_STROKE_HL 		= Graphics.getRGB(121,183,231);
var UI_FG_STROKE_SL 		= Graphics.getRGB(121,183,231);

var UI_TEXT 				= Graphics.getRGB(46,110,158);
var UI_TEXT_HL 				= Graphics.getRGB(29,89,135);
var UI_TEXT_SL 				= Graphics.getRGB(225,112,9);

/*
 * Fragment - inherits display object and represents a fragment which is in view
 *  
 * */
Fragment.prototype = new Shape();
Fragment.prototype.constructor = Fragment;
function Fragment(metadata)
{
	var self = this;
	Shape.call(this, new Graphics);
	this.meta = metadata;
	this.size = 0;//in degrees
	//this.rotation
	this.regX = 0; this.regY = 0;
	this.radius = 0;
	this.dragging = false;
	this.mouseOver = false;
	this.area = 'rm'; //or ccw or cw
	this.fill = UI_FG_FILL;
	this.stroke = UI_FG_STROKE;
	
	this.lower = function() {return this.rotation - this.size/2.0;};
	this.higher = function() {return this.rotation + this.size/2.0;};
	
	var anim = 0; //the number of animations currently going on, controls whether to redraw
	
	this.setAnim = function(){anim = anim+1;}
	this.clearAnim = function(){anim=anim-1; if(anim==0) self.redraw();}
	
	// -------------- Transitions
	this.setArea = function(a) 
	{
		if(self.area == a) return;
		self.area = a;
		switch(a)
		{
			case 'cw':
				self.x = 0; self.y = 0;
				self.radius = CW_R;
				break;
			case 'ccw':
				self.x = 0; self.y = 0;
				self.radius = CCW_R;
				break;			
		}
		if(self.dragging)
			self.radius = self.radius + DRAG_R;
	};
	this.setRadius = function(s) 
	{
		var target = {radius: Math.abs(s)};
		var tween = Tween.get(self)
			.call(setAnim).call(self.setAnim)
			.to(target, 250, Ease.quartOut)
			.call(self.clearAnim).call(clearAnim);
	};
	this.setRotation = function(a) {
		var r = bound_angle(a - this.rotation); //r is the distance and direction of shortest rotation
		
		var target = {rotation: this.rotation + r};
		var tween = Tween.get(self)
			.call(setAnim)
			.to(target, 250, Ease.quartOut)
			call(clearAnim).call('self.rotation = bound_angle(a);');
	};
	this.setSize = function(s)
	{
		var target = {size: Math.abs(s),};
		var tween = Tween.get(self)
			.call(setAnim).call(self.setAnim)
			.to(target, 250, Ease.quartOut)
			.call(self.clearAnim).call(clearAnim).call('self.rotation = bound_angle(a);');
	};
	this.setSizeBP = function(t) //t is total length of construct
	{
		this.setSize( 360.0 * (t / self.meta.length));
	}
	this.setDragging = function(d) //d = dragging
	{
		if(self.dragging == d) return;
		if(d)
		{
			self.setRadius(self.radius + DRAG_R);
		}
		else
		{
			self.setRadius(self.radius - DRAG_R);
		}
		self.dragging = d;
	}
	// -------------- end Transitions
	
	this.onMouseMove = function(x,y)
	{
		var ht = self.hitTest(x,y);
		if(ht && !self.mouseOver)
		{
			self.mouseOver = true;
			self.alpha = 0.8;
			_set_cursor('move');
		}
		if(!ht && self.mouseOver)
		{
			self.mouseOver = false;
			self.alpha = 1.0;
			_set_cursor('auto');
		}
	};
	
	this.onDrag = function(l) //l = local position
	{
		//console.log('  Fragment.onDrag ('+l.x+', '+l.y+') area: ' + self.area);
		switch(self.area)
		{
			case 'rm':
				self.x = l.x;
				self.y = l.y;
				break;
			case 'cw':
			case 'ccw':
				self.rotation = l.toRadial().t;
				break;
		}
		//console.log('  -- x,y = ('+self.x+', '+self.y+') -- r,a = ('+self.radius+', '+self.rotation+')');
		stage.update();
	}
	
	
	var _draw_frag = function(g, d) // d = direction E [1,-1], 1 -> CW, -1 -> CCW
	{
		g.moveToRA(self.radius-FRAG_W_2, -s2) //left inner
		 .lineToRA(self.radius, -s2 + d * FRAG_ARROW) //left center
		 .lineToRA(self.radius+FRAG_W_2, -s2) //left outer
		 .arc(0,0,self.radius+FRAG_W_2, -s2, s2, false) //outer edge
		 .lineTo(self.radius, s2 + d * FRAG_ARROW) //right center
		 .lineTo(self.radius-FRAG_W_2, s2) //right inner
		 .arc(0,0,self.radius-FRAG_W_2, s2, -s2, true); //inner edge
	}
	
	this.redraw = function()
	{
		console.log('Frag.redraw()');
		var g = self.graphics.clear();
		g.setStrokeStyle(1);
		g.beginFill(self.fill)
		 .beginStroke(self.stroke);
		var s2 = self.size / 2.0;
		switch(self.area)
		{
			case 'rm':
				g.rect(-FRAG_W * 3, -FRAG_W * .5, FRAG_W * 6, FRAG_W * .5);
				break;
			case 'cw':
				_draw_frag(g, +1);
				break;
			case 'ccw':
				_draw_frag(g, -1);
				break;
		}
	};
	this.tick = function()
	{
		if(anim > 0)
			self.redraw();
	}
	this.redraw();
};

/*
 * Item - a fragment as it appears in the library
 *  
 * */
Item.prototype = new Container();
Item.prototype.constructor = Item;
function Item(x,y,w,h,m) //x,y,w,h, metadata
{
	var self = this;
	Container.call(this);
	
	this.fill = UI_FG_FILL;
	this.stroke = UI_FG_STROKE;
	this.meta = m;
	this.x = x;
	this.y = y;
	
	var t = new Text(m.name, '500 12px Lucida Grande,Lucida Sans,Arial,sans-serif', UI_TEXT);
	t.textAlign = 'center';
	t.textBaseline = 'middle';
	t.x = w / 2.0;
	t.y = h / 2.0;
	var shape = new Shape(new Graphics);
	var r = new Rectangle(x,y,w,h);
	this.addChild(shape);
	this.addChild(t);
	
	this.onMouseMove = function(x,y) //return true if mouse is over me
	{
		if(r.containsXY(x,y))
		{
			//console.log('Mouse In: ' + self.meta.name);
			self.fill = UI_FG_FILL_HL;
			self.stroke = UI_FG_STROKE_HL;
			_set_cursor('move');
			self.mouseOver = true;
			_redraw();
			stage.update();
			return true;
		}
		else if(self.mouseOver)
		{
			//console.log('Mouse Out: ' + self.meta.name);
			self.fill = UI_FG_FILL;
			self.stroke = UI_FG_STROKE;
			_set_cursor('auto');
			self.mouseOver = false;
			_redraw();
			stage.update();
		}
		return false;
	};	
			
	var _redraw = function()
	{
		var g = shape.graphics;
		g.clear()
		 .beginFill(self.fill)
		 .beginStroke(self.stroke)
		 .drawRoundRect(0,0,w,h,4);		
	}
	
	_redraw();
};

/*
 * Library - the library
 *  
 * */
Library.prototype = new Container();
Library.prototype.constructor = Library;
function Library(x,y,w,h, slide)
{
	var self = this;
	Container.call(this);
	
	this.open = false;
	this.slide = slide;
	this.onSelect = null;
	this.x = x;
	this.y = y;
	this.width = w;
	this.height = h;
	
	var itemPanel;
	var itemHeight;
	
	//PUBLIC METHODS	
	this.Open = function() 
	{
		self.open = true;
		var target = {x: -self.slide,};
		var tween = Tween.get(slider)
			.call(setAnim)
			.to(target, 500, Ease.quartOut)
			.call(clearAnim);
	};
	this.Close = function()
	{
		self.open = false;
		var target = {x: 0,};
		var tween = Tween.get(slider)
			.call(setAnim)
			.to(target, 500, Ease.quartOut)
			.call(clearAnim);
	};
	this.addItem = function(m, update) //addItem(metadata, update = true)
	{
		var i = new Item(0,itemHeight,itemPanel.width,LIB_ITEM_HEIGHT,m);
		itemHeight = itemHeight + LIB_ITEM_HEIGHT + LIB_ITEM_SPACING;
		itemPanel.addChild(i);
		if(update == undefined) update = true;
		if(update)
			stage.update();
	};
	this.addItems = function(m) //addItems(Array[metadata])
	{
		for(i in m) self.addItem(m[i], false);
		stage.update();
	}
	this.onMouseMove = function(e)
	{
		//console.log('library.onMouseMove() ('+e.stageX+','+e.stageY+')');
		if(self.open)
		{
			if(!stage.mouseInBounds)
			{
				self.Close();
				return;
			}
			if(!rc.containsXY(e.stageX, e.stageY))
			{
				self.Close();
				return;
			}
			else
			{
				//are we moused over any items?
				var l = itemPanel.globalToLocal(e.stageX, e.stageY);
				for(var i = 0; i < itemPanel.getNumChildren(); i = i + 1)
				{
					itemPanel.getChildAt(i).onMouseMove(l.x,l.y);
				}
			}
		}
		if(!self.open)
		{
			if(ro.containsXY(e.stageX, e.stageY))
			{
				self.Open();
				return;
			}
		}
	};
	self.onMouseDown = function(e)
	{
		if(self.open)
		{
			for(var i = 0; i < itemPanel.getNumChildren(); i = i + 1)
			{
				if(itemPanel.getChildAt(i).mouseOver)
				{
					self.Close();
					self.onSelect(e, itemPanel.getChildAt(i).meta);
					return true;
				}
			}
		}
		return false;
	}
	//PRIVATE METHODS
	var _build_handle = function(handle_w, panel_w, h)
	{
		var g = new Graphics();
		//Draw the panel
		g.setStrokeStyle(1);
		g.beginFill(Graphics.getRGB(0,0,0,0.8));
		g.drawRect(handle_w, 8, panel_w, h - 16);
		
		//Draw the Handle
		g.beginFill(UI_BG_FILL);
		g.beginStroke(UI_BG_STROKE);
		g.drawRoundRectComplex( 0, 0, LIB_HANDLE_WIDTH, bounds.height, 4, 0, 0, 4);
		
		var s = new Shape(g);	
		
		//add the text
		var t = new Text('Library', '700 13.33px Lucida Grande,Lucida Sans,Arial,sans-serif', UI_FG_FILL);
		t.textBaseline = 'middle';
		t.rotation = -90;
		t.x = LIB_HANDLE_WIDTH / 2.0;
		t.y = bounds.height * 0.2; 
		
		var c = new Container();
		c.addChild(s);
		c.addChild(t);
		
		return c;
	};
		
	var _build_stopper = function(x,y,w,h)
	{
		var g = new Graphics();
		g.setStrokeStyle(1);
		g.beginFill(UI_BG_FILL);
		g.beginStroke(UI_BG_STROKE);
		g.drawRoundRectComplex( 0, 0, w, h, 0, 4, 4, 0);
		
		var s = new Shape(g);
		s.x = x;
		s.y = y;
		return s;
	};
	var _build_item_panel = function(x,y,w,h)
	{
		var r = new Container();
		r.x = x;
		r.y = y;
		r.height = h;
		r.width = w;
		return r;
	}
	
	
	//initiation object hirachy
	var hnd = _build_handle(LIB_HANDLE_WIDTH, self.slide, self.height);
	var stop = _build_stopper(LIB_HANDLE_WIDTH, 0, LIB_STOPPER_WIDTH, self.height);
	itemPanel = _build_item_panel(LIB_HANDLE_WIDTH + 4, 10, self.slide - 16, self.height - 10);
	itemHeight = LIB_ITEM_SPACING;
	
	var ro = new Rectangle(x,y,w,h);
	var rc = new Rectangle(x-this.slide,y,w+this.slide,h);
	
	var slider = new Container();
	slider.x = 0;
	slider.y = 0;
	slider.height = this.height;
	slider.width = this.slide;
	slider.addChild(hnd);
	slider.addChild(itemPanel);
			
	this.addChild(slider);
	this.addChild(stop);
	
}; //end Library def

/*
 * Designer - the main construct designer
 *  
 * */
Designer.prototype = new Container();
Designer.prototype.constructor = Designer;
function Designer(x,y,w,h) //x,y,w,h, metadata
{
	var self = this;
	Container.call(this);
	
	// ------------ Layout variables
	this.x = x; this.y = y; this.width = w; this.height = h;
	this.c = new Point((x+w) / 2.0, (y+h) / 2.0);
	CW_R = Math.min(h,w) * 0.5;
	CCW_R = Math.min(h,w) * 0.2;
	var ccw_limit = CCW_R + 1.5*FRAG_W;
	var cw_limit = CW_R + 1.5*FRAG_W;
	
	// --- display 
	var dzfill = Graphics.getRGB(60,60,60);
	var dz = {'rm':null, 'cw':null, 'ccw':null,}; //Dropzones
	var fr; //fragments
	var length = 0;
		
	// other vars
	var drag = null; //the element being dragged, null if none
	var active_dz = '';
	
	// Public methods
	
	this.initialFragments = function(m) //m = metadata from server
	{};
	this.addFragment = function(e, m) //m = metadata from library, e = mouseEvent
	{
		var f = new Fragment(m);
		f.size = 360 * (m.length / (m.length * length));
		f.dragging = true;
		drag = f;
		_set_active_dz('rm');
		e.onMouseMove = self.onFragDrag;
		e.onMouseUp = self.onFragDrop;
		fr.addChild(f);
		stage.update();
	};
	this.onFragDrag = function(e)
	{
		//console.log('designer.onFragDrag('+e.stageX+', '+e.stageY+')');
		drag.onDrag(fr.globalToLocal(e.stageX, e.stageY));//tell the fragment that we moved
		var d = _get_dz_under(e.mouseX, e.mouseY); 
		if(d != active_dz) //have we changed area?
		{
			_set_active_dz(d);
			drag.setArea(d);
		}
	};
	this.onFragDrop = function(e)
	{
		console.log('onFragDrop: ' + active_dz);
		switch(active_dz)
		{
			case 'rm':
				fr.removeChild(drag);
				drag = null;
				break;
			default:
				drag.setDragging(false);
				drag = null;
				break;
		}
		_set_cursor('auto');
		e.onMouseMove = null;
		e.onMouseUp = null;
		dz[active_dz].fadeOut();
		stage.update();
	};
	this.onMouseMove = function(e)
	{
		if(drag == null)
		{
			var l = fr.globalToLocal(e.stageX, e.stageY); 
			for(var f = 0; f < fr.getNumChildren(); f = f+1)
			{
				self.fr.getChildAt(f).onMouseMove(l.x,l.y);
			}
		}
	};
	this.onMouseDown = function(e)
	{};
	
	// Private methods
	var _init = function() 
	{
		console.log('designer._init() ' + ccw_limit + ' < ' + cw_limit);
		var d = new Container();
		dz.rm = new Shape(new Graphics());
		dz.cw = new Shape(new Graphics());
		dz.ccw = new Shape(new Graphics());
		d.addChild(dz.rm);
		d.addChild(dz.cw);
		d.addChild(dz.ccw);
		
		fr = new Container();
		fr.x = self.c.x;
		fr.y = self.c.y;
		
		self.addChild(d);
		self.addChild(fr);
		
		_draw_dz();
	};
	var _draw_dz = function() //draw all three drop zones
	{
		var c = self.c;
		//Remove zone
		dz.rm.graphics.clear()
		 .beginFill(dzfill)
		 .drawRoundRect(x,y,w,h,4)
		 .beginFill(Graphics.getRGB(255,255,255))
		 .drawCircle(c.x,c.y,cw_limit);
		dz.rm.alpha = 0.0;
		
		//clockwise zone
		dz.cw.graphics.clear()
		 .beginFill(dzfill)
		 .drawCircle(c.x,c.y,cw_limit)
		 .beginFill(Graphics.getRGB(255,255,255))
		 .drawCircle(c.x,c.y,ccw_limit);
		dz.cw.alpha = 0.0;
		
		//counter-clockwise zone
		dz.ccw.graphics.clear()
		 .beginFill(dzfill)
		 .drawCircle(c.x,c.y,ccw_limit);
		dz.ccw.alpha = 0.0;
	}
	var _set_active_dz = function(d) //set dropzone d ['rm', 'cw', 'ccw'] as active
	{
		if(active_dz == d) return;
		console.log('_set_active_dz("'+ d + '")');
		try
		{
			dz[active_dz].fadeOut(250);
		}
		catch(e) {};
		dz[d].fadeIn(250);
		if(active_dz == 'rm')
		{
			_set_length(length + drag.meta.length);
		}
		if(d == 'rm')
		{
			_set_length(length - drag.meta.length);
		}
		active_dz = d;
	};
	var _set_length = function(l)
	{
		if(l < 0) return;
		length = l;
		if(l>0)
		{
			for(var i = 0; i < fr.getNumChildren(); i = i +1)
			{
				fr.getChildAt(i).setLengthBP(l);
			}
		}
	};
	var _get_dz_under = function(x,y) //return dropzone under x,y as string
	{
		var d = self.c.distXY(x,y);
		if(d < self.ccw_limit)
			return 'ccw';
		if(d < self.cw_limit)
			return 'cw';
		else
			return 'rm';
	};
	
	_init();
};



var library;
var designer;

var $c;
//stage is only updated if an animation is playing
var anim = 0;
var setAnim = function() {anim=anim+1;};
var clearAnim = function() {anim=anim-1; if(anim == 0) stage.update();};

var init_designer = function(){	
	$c = $('#cdesigner');
	$c.prop('width', $c.width());
	$c.prop('height', $c.height());
	canvas = document.getElementById('cdesigner');
	stage = new Stage(canvas);
	//stage.enableMouseOver(25);
	
	_calc_size();
	
	
	library = new Library(LIB_X,LIB_Y,LIB_W,LIB_H,LIB_SLIDE);
	designer = new Designer(0,0,LIB_X - 2, bounds.height);
	library.onSelect = function(e, m)
	{
		console.log('Selected: ' + m.name);
		designer.addFragment(e,m);
	}
	stage.addChild(designer);
	stage.addChild(library);
	
	stage.onMouseMove = _on_mouse_move;
	stage.onPress = _on_press;
	Ticker.setFPS(25);
	Ticker.addListener(this);
	
	list_fragments(show_items);
	stage.update();
};
var tick = function()
{
	if(anim>0)
		stage.update();
};
var _calc_size = function(){
	bounds.width = canvas.width;
	bounds.height = canvas.height;
	LIB_X = Math.floor(bounds.width - LIB_HANDLE_WIDTH - LIB_STOPPER_WIDTH);
	LIB_Y = 0;
	LIB_SLIDE = Math.floor((1-0.618 ) * bounds.width)
	LIB_W = LIB_HANDLE_WIDTH + LIB_STOPPER_WIDTH;
	LIB_H = bounds.height;
};

var show_items = function(metadata)
{
	library.addItems(metadata);
	stage.update();
};

// Mouse things --------------------------
var _on_mouse_move = function(e)
{
	library.onMouseMove(e);
	if(!library.open) designer.onMouseMove(e);
};
var _on_press = function(e)
{
	if(library.open)
		library.onMouseDown(e);
	else
		designer.onMouseDown(e);
}

//Mouse Utils
var _is_over = function(ob, x,y)
{
	var l = ob.globalToLocal(x,y);
	return ob.hitTest(l.x,l.y);
}
var rp = Rectangle.prototype;
rp.containsPoint = function(p)
{
	return this.containsXY(p.x,p.y);
};
rp.containsXY = function(x,y)
{
	return (x > this.x && 
			x < (this.x + this.width) &&
			y > this.y &&
			y < (this.y + this.height));
};
var pp = Point.prototype;
pp.dist = function(p) //distance to another point
{
	return this.distXY(p.x,p.y);
};
pp.distXY = function(x,y)
{
	var a = x-this.x; var b = y-this.y;
	return Math.sqrt(a*a + b*b);
}
pp.angle = function(p) //degrees to another point from +ve X axis, CCW is +ve, limit +-180
{
	return this.angleXY(p.x,p.y);
}
var _R2D = 180.0 / Math.PI;
pp.angleXY = function(x,y)
{
	return Math.atan2(p.y-this.y,p.x-this.x) * _R2D;
}
pp.toRadial = function()
{
	return xy2ra(this.x,this.y);
}
var gp = Graphics.prototype;
gp.moveToRA = function(r,a)
{
	var p = ra2xy(r,a);
	return this.moveTo(p.x,p.y);
}
gp.lineToRA = function(r,a)
{
	var p = ra2xy(r,a);
	return this.lineTo(p.x,p.y);
}
var xy2ra = function(x,y)
{
	return {r: Math.sqrt(x*x+y*y), a: Math.atan2(y,x) * _R2D};
}
var ra2xy = function(r,a)
{
	var t = a / _R2D;
	return {x: r*Math.cos(t),y: r*Math.sin(t),};
}

var shortest_angle = function(a,b) //return the shortest distance between the two
{
	return Math.min(
					Math.abs(a-b), 
					Math.abs(360+b - a));
}
var bound_angle = function(a) //make a within [-180, 180]
{
	if(a > 180) return a - 360;
	if(a < -180) return a + 360;
	return a;
};

var dop = DisplayObject.prototype;
dop.fadeIn = function(t) //fade in in time t
{
	console.log('fadeIn');
	this.fade(t, 1.0);
};
dop.fadeOut = function(t) //fade in in time t
{
	console.log('fadeOut');
	this.fade(t, 0.0);
};
dop.fade = function(t, a) //fade to alpha a in in time t
{
	var target = {alpha: a};
	var tween = Tween.get(this)
		.call(setAnim)
		.to(target, t, Ease.quartOut)
		.call(clearAnim);
};

//eg auto, move, wait
var _set_cursor = function(type)
{
	if(type == undefined) type = 'auto';
	$c.css('cursor', type);
}

// prevent text cursor when dragging
window.addEventListener("dragstart", function(fEventObject){ CancelEvent(fEventObject); } );
window.addEventListener("selectstart", function(fEventObject){ CancelEvent(fEventObject); } );

function CancelEvent(fEventObject) 
{
   if (fEventObject.preventDefault) fEventObject.preventDefault();
   if (fEventObject.cancel != null) fEventObject.cancel = true;
}
