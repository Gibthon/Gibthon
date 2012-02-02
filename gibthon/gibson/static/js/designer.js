var canvas;
var stage;
var meta = null;

//Declare layout values
//constants
var LIB_HANDLE_WIDTH = 20;
var LIB_STOPPER_WIDTH = 6;
var LIB_ITEM_HEIGHT = 20;
var LIB_ITEM_SPACING = 4;

//variable
var bounds = new Rectangle();
var LIB_X = 0; var LIB_Y = 0; var LIB_W = 0; var LIB_H; var LIB_SLIDE = 0;

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

var BLACK = Graphics.getRGB(0,0,0);

var UI_FONT_LG = '700 24px Lucida Grande,Lucida Sans,Arial,sans-serif';
var UI_FONT_M = '500 12px Lucida Grande,Lucida Sans,Arial,sans-serif';
var UI_FONT_SM = '300 10px Lucida Grande,Lucida Sans,Arial,sans-serif';

var get_next_color = function()
{
	if(get_next_color.i == undefined)
		get_next_color.i = Math.floor(Math.random() * 256);
	get_next_color.i = (75 + get_next_color.i) % 255;
	var r = Graphics.getHSL(get_next_color.i,80,40);
	return r;
}

/*Areas:
 * 
 * Fragments can be in one of three areas as specified:
 * 
 * */
var CW = 0;
var CCW = 1;
var RM = 2;
//convert area to human readable string
var a2s = function(a) {switch(a){case CW: return 'cw'; case CCW: return 'ccw'; case RM: return 'rm'; default: return 'unknown';}}

/*
 * Item - a fragment as it appears in the library
 * 
 * CONDEMNED! -- library to become a jquery widget
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
 * CONDEMNED! -- library to become a jquery widget
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
 * Label - show a display fragment's label
 *  
 * */
Label.prototype = new Container();
Label.prototype.constructor = Label;
function Label(text, radius, font, color)
{
	var self = this;
	Container.call(this);
	
	if(font == undefined) font = UI_FONT_SM;
	if(color == undefined) color = BLACK;
	
	var t = text;
	var r = radius;
	var f = font;
	var col = color;
	var s = 1; //controls which way up the text should display +1 => bottom inwards, -1 => top inwards
	
	//holds the text
	var c = new Container();
	self.addChild(c);
	
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++ Public methods
	
	/*
	 * update Rotation - should the text face inwards or outwards
	 * 
	 * a: the absolute angle relative to +ve x axis
	 * */
	this.updateRotation = function(a, update)
	{
		var x = Math.sin(a / _R2D);
		var _s = s;
		if(x > 0)
			_s = -1;
		if(x < 0)
			_s = 1;
		if(s != _s)
		{
			s = _s;
			reposition();
		}
		
		_u(update);
	};
	
	this.setRadius = function (_r, update)
	{
		if(_r != r)
		{
			r = _r;
			reposition();
			_u(update);
		}
	};
	
	this.setText = function(_t, update)
	{
		if(_t != t)
		{
			r = _r;
			redraw();
			_u(update);
		}
	};
	
	this.getRadius = function() {return r;};
	this.getText = function() {return t;};
	
	this.show = function(update) 
	{
		self.visible = true;
		_u(update);
	};
	this.hide = function(update) 
	{
		self.visible = false;
		_u(update);
	};
	
	// ------------------------------------------------------ Private methods
	
	var reposition = function()
	{
		var a = 0;
		//for each letter
		for(var i = 0; i < c.getNumChildren(); i = i + 1)
		{
			var l = c.getChildAt(i);
			l.regY = s * r;
			l.x=0;
			l.rotation = a;
			a = a + (s * (l.getMeasuredWidth() + 0.5)/ r) * _R2D;
		}
		c.rotation = bound_angle(- a / 2.0 + s * 90);
		self.uncache();
		var _r = r + 10;
		self.cache(-_r,-_r, 2 * _r, 2 * _r);
	};
	
	var redraw = function()
	{
		//clear any previous text
		c.removeAllChildren();
		
		var a = 0;
		//for each letter
		for(var i = 0; i < t.length; i = i + 1)
		{
			var l = new Text(t[i], f, col);
			l.textBaseline = 'middle';
			c.addChild(l);
		}
		reposition();
	};
	
	var _u = function(u)
	{
		if(u == undefined) u = true;
		if(u) stage.update();
	}
	
	//intiate :
	redraw();
	
}; //end Label


/*	F - global settings for DisplayFragments
 * 
 * 
 * */
var F = new Object();
F.radii = [0,0,0]; //radii at which to render for each area
F.label = [0,0,0]; //radii at which to render the labels
F.width = 25; // width of fragment
F.arrow = 5; //how many pixels the arrow should be when in the construct
F.delta = 5; //how much of a gap to leave when dragging

/*
 * DisplayFragment - inherits display object and represents a fragment which is in view in the Designer
 *  
 * */
DisplayFragment.prototype = new Container();
DisplayFragment.prototype.constructor = Fragment;
//DisplayFragment(Fragment f, ConstructFragment cf = undefined, Area _area = RM)
function DisplayFragment(f, cf)
{
	var self = this;
	Shape.call(this, new Graphics);
	
	/*
	 * this.f - stores the fragment related to this DisplayFragment - go here for metadata
	 * */
	this.f = f;
	
	/*
	 * this.cf - stores the constructFragment for this object.
	 * 
	 * */
	this.cf = cf;
	
	// Display things
	this.area = RM;
	if(cf.id != undefined)
	{
		if(cf.strand == 1) this.area = CW;
		else this.area = CCW;
	}
	
	this.angle = 0; //in degrees
	
	this.radius = F.radii[this.area];
	
	this.mouseOver = false;
	
	this.color = get_next_color();//UI_BG_FILL;
	this.stroke = Graphics.getRGB(255,255,255);
	this.drag = false;
	
	this.regX = 0; this.regY = 0;
	
	this.start = function() {return this.rotation;};
	this.middle = function() {return this.rotation + self.angle/2.0;};
	this.end = function() {return this.rotation + self.angle;};
	
	var anim = 0; //the number of animations currently running, controls whether to redraw on tick
	
	this.setAnim = function()
	{
		anim = anim+1;
		if(this.tick == undefined)
			this.tick = self.redraw;
	}
	this.clearAnim = function()
	{
		anim=anim-1;
		if(anim==0)
		{
			self.tick = undefined;
			self.redraw();
		}
	}
	
	/*
	 * my label
	 * */
	var l = new Label(self.f.name, F.label[self.area]);
	if(this.area == RM)
		l.hide();
	self.addChild(l);
	
	var g = new Graphics();
	var s = new Shape(g);
	self.addChild(s);
	
	this.setAngle = function(a)
	{
		self.angle = a;
		l.rotation = a / 2.0;
		l.updateRotation(self.middle());
	}
	this.getAngle = function() {return angle;}
	
	//set the fragment's area -- not animated
	this.setArea = function(a) 
	{
		if(self.area == a) return;
		var a_ = self.area;
		self.area = a;
		if(a == CW || a == CCW)
		{
			self.x = 0; self.y = 0;
			l.setRadius(F.label[a]);
			l.show();
		}
		else(a == RM)
			l.hide();
		
		this.radius = F.radii[a];
		if(this.drag)
		{
			l.hide();
			if(a == CW) this.radius = this.radius + F.delta;
			if(a == CCW) this.radius = this.radius - F.delta;
		}
		self.redraw();
	};
	
	this.setDrag = function(d, a, cb)//a is angle for end drag
	{
		if(d == this.drag) return;
		if(d) _start_drag();
		else _end_drag(a, cb);
	};
	
	var _start_drag = function(cb)
	{
		if(self.area == CW) self.radius = F.radii[CW] + F.delta;
		if(self.area == CCW) self.radius = F.radii[CCW] - F.delta;
		self.drag = true;
		self.redraw();
		l.hide();
		if(cb != undefined) cb();
	};
	
	var _end_drag = function(a)
	{
		var t = {'radius': F.radii[self.area],'rotation':a,};
		var cb = function() {l.show();};
		self.animate(t, cb);
		self.drag = false;
	};
	
	this.animate = function(t, cb) //animate angle, rotation or radius
	{
		if(t.rotation != undefined)
		{
			var r = bound_angle(t.rotation - this.rotation); //r is the distance and direction of shortest rotation
			t.rotation = this.rotation + r;
		}
		//make sure that the rotation is bound at the end of the animation, and make sure the label keeps up
		var done = function() {self.rotation = bound_angle(self.rotation); self.setAngle( self.angle );};
		
		//r: whether we need to redraw at each frame
		var r = (t.angle != undefined) || (t.radius != undefined);
		var tween = Tween.get(self).call(setAnim);
		if(r) tween.call(self.setAnim);
		tween.to(t, 250, Ease.quartOut)
		if(r) tween.call(self.clearAnim)
		tween.call(done).call(clearAnim);
		if(cb != undefined) tween.call(cb);
	};

	this.onMouseMove = function(r)
	{
		var ht = rcontains(r);
		if(ht && !self.mouseOver)
		{
			self.mouseOver = true;
			self.alpha = 0.8;
			_set_cursor('move');
			stage.update();
		}
		if(!ht && self.mouseOver)
		{
			self.mouseOver = false;
			self.alpha = 1.0;
			_set_cursor('auto');
			stage.update();
		}
		return self.mouseOver;
	};
	
	var rcontains = function(r) //does the fragment contain the radial point r
	{
		if( (r.r > self.radius + F.width / 2.0) || (r.r < self.radius - F.width / 2.0))
		{
			return false;
		}
		return a_contains(self.rotation, self.end(), r.a);
	};
	
	var _draw_frag = function(d) // d = direction E [1,-1], 1 -> CW, -1 -> CCW
	{
		var w2 = F.width / 2.0;
		var l = self.angle / _R2D;
		var a = F.arrow / self.radius;
		g.moveToRA(self.radius - w2, 0) //left inner
		 .lineToRA(self.radius, d * a) //left center
		 .lineToRA(self.radius + w2, 0) //left outer
		 .arc(0,0,self.radius + w2, 0, l, false) //outer edge
		 .lineToRA(self.radius +w2, l ) //right outer
		 .lineToRA(self.radius, l + d * a) //right center
		 .lineToRA(self.radius-w2, l) //right inner
		 .arc(0,0,self.radius - w2, l, 0, true) //inner edge
		 .closePath();
	};
	
	this.redraw = function()
	{
		g.clear()
		 .setStrokeStyle(1)
		 .beginFill(self.color)
		 .beginStroke(self.stroke);

		switch(self.area)
		{
			case RM:
				g.rect(-F.width * 3, -F.width * .5, F.width * 6, F.width * .5);
				break;
			case CW:
				_draw_frag(+1);
				break;
			case CCW:
				_draw_frag(-1);
				break;
		}
		
	};
	
	this.tick = undefined;
	
	this.redraw();
};

/*
 * 
 * 
 * Designer - the main construct designer
 *  
 * */
Designer.prototype = new Container();
Designer.prototype.constructor = Designer;
function Designer(x,y,w,h,cid) //x,y,w,h, metadata
{
	var self = this;
	Container.call(this);
	
/*
 * Designer Variables
 * */
 
	/*
	 * Designer - construct properties
	 * */
	this.name = '';
	this.len = 0;
	this.cid = cid;
	
	/*
	 * Designer layout and sizing variables
	 * */
	var c = new Point(w*0.5,h*0.5);
	var radius = Math.min(w,h) * 0.35; //the base radius of the plasmid
	
	F.width = radius * 0.1;
	F.delta = radius * 0.15;
	F.radii[CW] = radius + F.delta;
	F.radii[CCW]= radius - F.delta;
	F.radii[RM] = 1000;
	
	F.label[CW] = radius + 1.5 * F.delta;
	F.label[CCW]= radius - 1.5 * F.delta;
	F.radii[RM] = 1000;
	 
	var title_font = UI_FONT_LG;
	var len_font = UI_FONT_M;
	var drag = false; //is a fragment being dragged
	var selected = null; //which fragment is selected
	
	/*
	 * Designer Display Objects
	 * */
	var name_t, len_t;
	var fc = new Container(); //fragment container
	
/*
 * Designer Functions
 * */
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++PUBLIC++++++++++
	this.setName = function(n) //set the plasmid's name
	{
		self.name = n;
		name_t.text = n;
		stage.update();
	};
	
	this.addFragment = function(e,f) //add a fragment
	{
		var f = new DisplayFragment(f,new ConstructFragment({}, f));
		var p = self.globalToLocal(e.stageX,e.stageY);
		
		f.x = p.x;
		f.y = p.y;
		f.drag = true;
		drag = true;

		e.onMouseMove = function(e) {self.dragFragment(e, f);};
		e.onMouseUp   = function(e) {self.dropFragment(e, f);};
		
		self.addChild(f);
	};
	
	this.removeFragment = function(df)
	{
		var g = df.cf.order;
		
		fc.removeFragAt(g);
		_set_len(self.len - df.cf.length());
		
		cd_rm_fragment(function() {}, self.cid, df.cf.id);
		
		g = b(g-1);
		console.log('update_layout because a fragment was removed');
		_update_layout(g);
	}
	
	this.dragFragment = function(e, f) //fragment moved
	{
		var p = fc.globalToLocal(e.stageX,e.stageY);
		var rp = p.toRadial();
		
		var a_ = _get_area(rp);
		//if the fragment's area changed?
		if(a_ != f.area)
		{
			//if we're joining the construct
			if(f.area == RM)
			{
				_join(f, rp.a);
			}
			if(a_ == RM)
			{
				_leave(f);
			}
			f.setArea(a_);
		}
		
		if(f.area == RM)
		{
			f.x = p.x + c.x;
			f.y = p.y + c.y;
		}
		else
		{
			f.rotation = rp.a + f._mo - f.angle / 2.0;
			
			//Check if we crossed any fragments
			if(fc.getNumChildren() > 2 && !a_contains(f._p, f._n, f.middle()) && (anim <= 0))
			{
				console.log('Crossed: a_contains('+f._p+', '+f._n+', '+f.middle()+') = false');
				var i = fc.getChildIndex(f);
				var i_ = i;
				var u = i;
				// if we're nearest the previous fragment
				if(Math.abs(shortest_angle(rp.a,f._p)) < Math.abs(shortest_angle(rp.a,f._n)))
				{
					i_ = b(i - 1);
					u = b(i + 1);
				}
				else
				{
					i_ = b(i + 1);
					u = b(i - 1);
				}
					
				//switch the fragments
				fc.swap(i,i_);
				console.log('update layout because fragment '+i+' swapped with '+i_);
				_update_layout(u);
			}
		}
		stage.update();		
	}
	
	this.dropFragment = function(e, f) //fragment dropped
	{
		_set_cursor();
		
		drag = false;
		if(f.area == RM)
		{
			self.removeChild(f);
			stage.update();
			
			if(f.cf.id != undefined)
			{
				cd_rm_fragment(function() {}, self.cid, f.cf.id);
			}
			return;
		}
		var o = fc.getChildIndex(f);
		f.setDrag(false, fc.getChildAt(b(o-1)).end());
		
		//let the server know what's going on
		if(f.cf.id == undefined)
		{
			cd_add_fragment(function(cf) {f.cf = new ConstructFragment(cf);},
							self.cid,
							f.f.fid,
							f.f.order);
			return;
		}
		
		f.cf.order = o;
		cf = new Array();
		d = new Array();
		for(var i = 0; i < fc.getNumChildren(); i = i + 1)
		{
			var f = fc.getChildAt(i);
			cf.push(f.cf.id);
			if(f.area == CW)
			{
				d.push(1);
			}
			else
			{
				d.push(-1);
			}
		}
		cd_reorder_fragments(function() {}, self.cid,{'cfid': cf, 'direction': d,});
	}
	
	this.onMouseMove = function(e)
	{
		if(!drag)
		{
			var r = fc.globalToLocal(e.stageX,e.stageY).toRadial();
			selected = null;
			for(var i = 0; i < fc.getNumChildren(); i = i + 1)
			{
				if(fc.getChildAt(i).onMouseMove(r))
					selected = i;
			}
		}
	};
	
	this.onMouseDown = function(e)
	{
		if(selected != null)
		{
			var s = selected;
			var t = setTimeout(function() {self.onMouseDrag(e,s);}, 250);
			e.onMouseUp = function() {clearTimeout(t); self.onMouseClick(e,s);};
		}
	};
	
	//a click to the fragment at s
	this.onMouseClick = function(e, s)
	{
		var df = fc.getChildAt(s);
		var loc = ra2xy(df.radius, df.middle());
		var g = fc.localToGlobal(loc.x, loc.y);
		show_info(df, g.x,g.y);
	}
	
	//begins dragging the fragment at s
	this.onMouseDrag = function(e, s)
	{
		e.onMouseUp = function() {};
		var f = fc.getChildAt(s);
		var r = fc.globalToLocal(e.stageX,e.stageY).toRadial();
		var t = {};
		
		drag = true;

		f._p = fc.getChildAt(b(s - 1)).middle();
		f._n = fc.getChildAt(b(s + 1)).middle();
		f._mo = bound_angle(f.middle() - r.a);
		
		f.setDrag(true);
		
		e.onMouseMove = function(e) {self.dragFragment(e,f);};
		e.onMouseUp = function(e) {self.dropFragment(e,f);};
	}
	
	this.construct2Local = function(r,t,target) //translate to global coords
	{
		var p = ra2xy(r, t);
		return fc.localToLocal(p.x,p.y,target);
	}; 
	
	// -----------------------------------------------------------PRIVATE---------
	//add DisplayFragment f to the construct at angle a and remove from self
	var _join = function(df, a) 
	{
		var i = _closest_gap(a);
		df._mo = 0; //mouse offset
		self.removeChild(df);
		fc.addFragAt(df, i);
		df.setAngle(360 * (df.cf.length() / (self.len + df.cf.length())));
				
		_set_len(self.len + df.cf.length());
		
		console.log('update layout because a fragmen joined');
		_update_layout(b(i-1));
	}
	
	var _leave = function(df) //remove fragment f to the construct and add to self
	{
		var g = df.cf.order;
		df._mo = undefined;
		fc.removeFragAt(g);
		df.rotation = 0;
		self.addChild(df);
		_set_len(self.len - df.cf.length());
		
		g = b(g-1);
		
		console.log('update layout because a fragment left');
		_update_layout(g);
	}
	
	var b = function(j)//make sure j is within fc array bounds
	{
		if(j < 0) j = j + fc.getNumChildren();
		return j % fc.getNumChildren();
	};
	
	var _closest_gap = function(a) //return the index of the fragment after closest gap to the angle given
	{
		var n = fc.getNumChildren();
		if(n == 0) return 0;
		a = bound_angle(a);
		
		/* -- I can do better...
		// 'static' var i is the initial guess
		if ( typeof _closest_gap.i == 'undefined' ) {
			_closest_gap.i = 0;
		}
		* */
		
		//var i = b(_closest_gap.i);
		
		//find the closest
		var d = 360;
		var c = -1;
		
		for(var i = 0; i < n; i = i + 1)
		{
			var d_ = Math.abs(shortest_angle(a, fc.getChildAt(i).rotation));
			if(d_ < d)
			{
				d = d_;
				c = i;
			}
		}

		return c;
		
	}
	
	var _update_layout = function(s, a) //layout the fragments tip to tail. start with index s at angle a
	{
		var n = fc.getNumChildren();
		if(n == 0) 
		{
			stage.update();
			return;
		}
		if(s==undefined) s = 0;

		s = b(s);

		if(a==undefined) a = fc.getChildAt(s).rotation;
		a = bound_angle(a);
		
		//calculate the effective length
		var min = self.len * 5.0 / 360.0;
		var el = 0;
		var len = new Array();
		for(var i = 0; i < n; i = i + 1)
		{
			var l = fc.getChildAt(i).cf.length();
			if(l > min)
				el = el + l;
			else
				el = el + min;
			len.push(l);
		}
		
		
		var c = a;
		var d = -1;
		var targets = new Array();
		for(var i = 0; i < n; i = i + 1)
		{
			var j = b(i + s);
			var f = fc.getChildAt(j);
			
			var l = 0;
			if(len[j] > min)
				l = 360.0 * (len[j] / el);
			else
				l = 360.0 * (min / el);
			
			if(!f.drag)
				targets.push({rotation:c,angle:l,});
			else
			{
				d = i;
				targets.push({});
			}
			c = bound_angle(c + l);
		}
		
		if(d >= 0)
		{
			var f = fc.getChildAt(b(d+s));
			var nxt = b(d+1);
			var pre = b(d-1);
			f._n = bound_angle(targets[nxt].rotation + targets[nxt].angle / 2.0);
			f._p = bound_angle(targets[pre].rotation + targets[pre].angle / 2.0);
		}
		
		for(var i = 0; i < n; i = i+1)
		{
			var j = b(i + s);
			fc.getChildAt(j).animate(targets[i]);
		}
	}
	
	var _get_area = function(r) //return the area based on radial point r
	{
		if(r.r < radius) return CCW;
		if(r.r < radius + F.delta + 3 * F.width) return CW;
		return RM;
	}
	
	var _set_len = function(l) //set the length
	{
		self.len = Math.abs(parseInt(l));
		len_t.text = self.len + ' bp';
	}
	
	var _dbg_str = function()
	{
		var s = 'Designer: ';
		for(var i = 0; i < fc.getNumChildren(); i = i + 1)
		{
			var c = fc.getChildAt(i);
			s = s +'\n  ' + i + ') ' + c.f.name + ' ['+c.cf.length()+'] ';
		}
		return s;
	}
	
	var _init = function(m)
	{
		self.name = m.name;
		self.len = 0;
		
	//draw name and length
		//setup the text
		name_t = new Text(self.name, title_font, BLACK);
		len_t = new Text(self.len + ' bp', len_font, BLACK);
		name_t.textAlign = 'center';
		len_t.textAlign = 'center';
		name_t.textBaseline = 'middle';
		len_t.textBaseline = 'middle';
		len_t.y = 25;
		name_t.maxWidth = 2*radius - 10; //CHANGE ME!
		len_t.maxWidth = 2*radius - 10; //CHANGE ME!
		
		//draw the plasmid
		var g = new Graphics();
		g.setStrokeStyle(1)
		 .beginStroke(BLACK)
		 .drawCircle(0,0,radius);
		var s = new Shape(g);
		
		//crate the container, set to center
		var con = new Container();
		con.x = c.x; 
		con.y = c.y;
			
		//add children
		con.addChild(s);
		con.addChild(len_t);
		con.addChild(name_t);
		
		fc.x = c.x;
		fc.y = c.y;
		
		//additional functions for fc
		fc.addFragAt = function(df, i) //add a DisplayFragment at position i
		{
			df.cf.order = i;
			this.addChildAt(df, i);
			for(var j = i + 1; j < this.getNumChildren(); j = j+1)
			{
				this.getChildAt(j).cf.order = j;
			}
		};
		fc.removeFragAt = function(i) //remove a DisplayFragment from position i
		{
			this.removeChildAt(i);
			for(var j = i; j < this.getNumChildren(); j = j+1)
			{
				this.getChildAt(j).cf.order = j;
			}
		};
		fc.removeFrag = function(df)
		{
			this.removeFragAt(this.getChildIndex(df));
		}
		fc.swap = function(i, j) //swap the order of the DisplayFragments at i and j
		{
			if(i == j) return;
			if(i > j)
			{
				var _i = j;
				j = i;
				i = _i;
			}
			var fj = this.getChildAt(j);
			var fi = this.getChildAt(i);
			fj.cf.order = i;
			fi.cf.order = j;
			this.sortChildren(function(a,b) {return a.cf.order - b.cf.order;});
		};
		
		//initiate the fragments
		for(var i in m.cfs)
		{
			var f = meta[m.cfs[i].fid];
			var cf = new ConstructFragment(m.cfs[i], f);
			var df = new DisplayFragment(f,cf);
			fc.addChild(df);
			self.len = self.len + cf.length();
		}
		
		
		self.addChild(con);
		self.addChild(fc);
		
		_set_len(self.len);
		console.log('initial update layout');
		_update_layout();
	};
	
	cd_get_info(_init, this.cid);
};

/*
 * 
 * 
 * 
 * 
 *  Initiation 
 * 
 *
 *
 *
 *
 *
 */


var library;
var designer;

var $c;
var $info;

//stage is only updated if an animation is playing
var anim = 0;
var setAnim = function() {anim=anim+1;};
var clearAnim = function() {anim=anim-1; if(anim == 0) stage.update();};

var init_designer = function(cid){
	$c = $('#cdesigner');
	init_info();
	var w = $c.width();
	var h = (7.0 / 16.0) * w;
	$c.height(h);
	$c.prop('width', w);
	$c.prop('height', h);
	
	
	canvas = document.getElementById('cdesigner');
	stage = new Stage(canvas);
	
	list_fragments(function(m) {continue_init(m,cid);});
	_calc_size();
	
	library = new Library(LIB_X,LIB_Y,LIB_W,LIB_H,LIB_SLIDE);
	library.onSelect = function(e, m)
	{
		designer.addFragment(e,m);
	}

	Ticker.setFPS(25);
	Ticker.addListener(this);
	
	stage.update();
};
var continue_init = function(m, cid)
{
	library.addItems(m);
	meta = {};
	for(var i in m)
	{
		meta[m[i].fid] = m[i];
	}
	designer = new Designer(0,0,LIB_X, bounds.height, cid);
	stage.addChild(designer);
	stage.addChild(library);
	
	stage.onMouseMove = _on_mouse_move;
	stage.onPress = _on_press;
	
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


// --------------------------------------Info window
var init_info = function()
{
	$info = $('.fragment-info');
	
	//setup the buttons
	$info.find('#fragment_remove')
		.button({label: 'Remove', icons: {primary:'ui-icon-trash',},})
		.click(function() {hide_info()});

	$info.find('#fragment_clip')
		.button({label: 'Clipping', icons: {primary:'ui-icon-scissors'}, disabled:true,})
		.click(function() {});
}

var show_info = function(df, x, y)
{
	//set the title and description
	$info.find('#fragment_name').text(df.f.name);
	$info.find('#fragment_desc').text(df.f.desc);
	
	//set position
	$info.css({position: 'absolute', zindex:100, left:x - 37, top:y - $info.height(),});
	$info.css({display:'block',});
	
	//set remove callback
	$info.find('#fragment_remove')
		.unbind('click')
		.click(function() {hide_info();designer.removeFragment(df);});
	
	//binding to click means it gets triggered immediately
	$c.mousedown(hide_info);
}

var hide_info = function()
{
	$info.css({display:'none',});
	$c.unbind('mousedown', hide_info);
}

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
	var p = ra2xy(r,a * _R2D);
	return this.moveTo(p.x,p.y);
}
gp.lineToRA = function(r,a)
{
	var p = ra2xy(r,a * _R2D);
	return this.lineTo(p.x,p.y);
}
gp.arcToRA = function(r, a1, a2)
{
	var p1 = ra2xy(r,a1);
	var p2 = ra2xy(r,a2);
	return this.arcTo(p1.x,p1.y,p2.x,p2.y,r);
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
var r2xy = function(r)
{
	var t = r.a / _R2D;
	return {x: r.r*Math.cos(t),y: r.r*Math.sin(t),};
}

var sign = function(a) {if(a==0) return 0; return (a / Math.abs(a));};
var shortest_angle = function(a,b) //return the shortest distance from a to b
{
	return bound_angle(b-a);
}
var bound_angle = function(a) //make a within [-180, 180]
{
	if(a > 180) return a - 360;
	if(a < -180) return a + 360;
	return a;
};
var a_contains = function(l,h,a) //is angle a in the segment defined by moving clockwise from l to h?
{
	if(isNaN(l) || isNaN(h) || isNaN(a)){
		console.log('a_contains('+l+','+h+','+a+')');
		return true;
	}
	var h_ = h; var a_ = a;
	if(h_ < l) h_ = h_ + 360;
	if(a_ < l) a_ = a_ + 360;
	return (a_ - l) < (h_ - l);
}
var dop = DisplayObject.prototype;
dop.fadeIn = function(t) //fade in in time t
{
	//console.log('fadeIn');
	this.fade(t, 1.0);
};
dop.fadeOut = function(t) //fade in in time t
{
	//console.log('fadeOut');
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
