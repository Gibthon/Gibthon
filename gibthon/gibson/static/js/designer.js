var canvas;
var stage;
var meta = null;

//Declare layout values
//constants
var LIB_HANDLE_WIDTH = 20;
var LIB_STOPPER_WIDTH = 6;
var LIB_ITEM_HEIGHT = 20;
var LIB_ITEM_SPACING = 4;

var FRAG_W = 25;
var FRAG_W_2 = FRAG_W / 2.0;
var FRAG_ARROW = 10.0;

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

var BLACK = Graphics.getRGB(0,0,0);

var UI_FONT_LG = '700 24px Lucida Grande,Lucida Sans,Arial,sans-serif';
var UI_FONT_M = '500 12px Lucida Grande,Lucida Sans,Arial,sans-serif';
var UI_FONT_SM = '300 10px Lucida Grande,Lucida Sans,Arial,sans-serif';

//Areas:
var CW = 0;
var CCW = 1;
var RM = 2;

var a2s = function(a) {switch(a){case CW: return 'cw'; case CCW: return 'ccw'; case RM: return 'rm'; default: return 'unknown';}}

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
 * Fragment - inherits display object and represents a fragment which is in view in the Designer
 *  
 * */
Fragment.prototype = new Shape();
Fragment.prototype.constructor = Fragment;
function Fragment(metadata, _area, w, r)
{
	var self = this;
	Shape.call(this, new Graphics);
	
	if(_area == undefined) _area = RM;
	if(w == undefined) w = 20;
	if(r == undefined) r = 450;
		
	this.meta = metadata;
	this.length = 0;//in degrees
	this.width = w;
	this.regX = 0; this.regY = 0;
	this.radius = r;
	this.mouseOver = false;
	this.area = _area; //or ccw or cw
	this.fill = UI_BG_FILL;
	this.stroke = Graphics.getRGB(255,255,255);
	this.drag = false;
	
	this.center = function() {return this.rotation + this.length/2.0;};
	this.end = function() {return this.rotation + this.length;};
	
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
			case CW:
				self.x = 0; self.y = 0;
				break;
			case CCW:
				self.x = 0; self.y = 0;
				break;			
		}
		self.redraw();
	};
	
	this.animate = function(t) //animate length, rotation, radius or total_length
	{
		var done = function() {};
		if(t.rotation != undefined)
		{
			var r = bound_angle(t.rotation - this.rotation); //r is the distance and direction of shortest rotation
			t.rotation = this.rotation + r;
			done = function() {self.rotation = bound_angle(t.rotation);};
		}

		if(t.total_length != undefined)
		{
			t.length = 360.0 * (self.meta.length / t.total_length);
			t.total_length = undefined;
		}
		var r = (t.length != undefined) || (t.radius != undefined);
		var tween = Tween.get(self).call(setAnim);
		if(r)
			tween.call(self.setAnim);
		tween.to(t, 250, Ease.quartOut)
		if(r)
			tween.call(self.clearAnim)
		tween.call(clearAnim).call(done);
	};

	// -------------- end Transitions
	
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
		if( (r.r > self.radius + self.width / 2.0) || (r.r < self.radius - self.width / 2.0))
		{
			return false;
		}
		return a_contains(self.rotation, self.end(), r.a);
	}	
	var _draw_frag = function(g, d) // d = direction E [1,-1], 1 -> CW, -1 -> CCW
	{
		var w2 = self.width / 2.0;
		var l = self.length / _R2D;
		var a = FRAG_ARROW / self.radius;
		g.moveToRA(self.radius - w2, 0) //left inner
		 .lineToRA(self.radius, d * a) //left center
		 .lineToRA(self.radius + w2, 0) //left outer
		 .arc(0,0,self.radius + w2, 0, l, false) //outer edge
		 .lineToRA(self.radius +w2, l ) //right outer
		 .lineToRA(self.radius, l + d * a) //right center
		 .lineToRA(self.radius-w2, l) //right inner
		 .arc(0,0,self.radius - w2, l, 0, true) //inner edge
		 .closePath();
	}
	
	this.redraw = function()
	{
		var g = self.graphics.clear();
		g.setStrokeStyle(1)
		 .beginFill(self.fill)
		 .beginStroke(self.stroke);

		switch(self.area)
		{
			case RM:
				g.rect(-self.width * 3, -self.width * .5, self.width * 6, self.width * .5);
				break;
			case CW:
				_draw_frag(g, +1);
				break;
			case CCW:
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
	var fragw = radius * 0.1;
	var delta = radius * 0.15;
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
	
	this.addFragment = function(e,m) //add a fragment
	{
		var f = new Fragment(m,RM,fragw);
		var p = self.globalToLocal(e.stageX,e.stageY);
		
		f.x = p.x;
		f.y = p.y;
		f.drag = true;
		drag = true;

		e.onMouseMove = function(e) {self.dragFragment(e, f);};
		e.onMouseUp   = function(e) {self.dropFragment(e, f);};
		
		self.addChild(f);
	};
	
	this.dragFragment = function(e, f) //fragment moved
	{
		var p = fc.globalToLocal(e.stageX,e.stageY);
		var rp = p.toRadial();
		
		var a_ = _get_area(rp);
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
			if(a_ == CW)
				f.radius = radius + delta + fragw;
			if(a_ == CCW)
				f.radius = radius - delta - fragw;
				
			f.setArea(a_);
		}
		
		if(f.area == RM)
		{
			f.x = p.x + c.x;
			f.y = p.y + c.y;
		}
		else
		{
			f.rotation = rp.a + f._mo - f.length / 2.0;
			
			//Check if we crossed any fragments
			if(fc.getNumChildren() > 2 && !a_contains(f._p, f._n, f.center()))
			{
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
				fc.removeChildAt(i);
				fc.addChildAt(f, i_);
				f._i = i_;
				_update_layout(u);
			}
		}
		stage.update();		
	}
	
	this.dropFragment = function(e, f) //fragment dropped
	{
		_set_cursor();
		f.drag = false;
		drag = false;
		if(f.area == RM)
		{
			self.removeChild(f);
			stage.update();
			console.log(_dbg_str());
			return;
		}
		var i = fc.getChildIndex(f);
		t = {rotation: fc.getChildAt(b(i-1)).end(),};
		if(f.area == CW)
		{
			t.radius = radius + delta;
		}
		if(f.area == CCW)
		{
			t.radius = radius - delta;
		}
		console.log(_dbg_str());
		
		f.animate(t);
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
			f = fc.getChildAt(selected);
			r = fc.globalToLocal(e.stageX,e.stageY).toRadial();
			t = {};
			f.drag = true;
			drag = true;
			f._p = fc.getChildAt(b(selected - 1)).center();
			f._n = fc.getChildAt(b(selected + 1)).center();
			f._mo = bound_angle(f.center() - r.a);
			
			if(f.area == CW)
				t.radius = radius + delta + fragw;
			if(f.area == CCW)
				t.radius = radius - delta - fragw;
				
			f.animate(t);
			
			e.onMouseMove = function(e) {self.dragFragment(e,f);};
			e.onMouseUp = function(e) {self.dropFragment(e,f);};
		}
	}
	
	// -----------------------------------------------------------PRIVATE---------
	var _join = function(f, a) //add fragment f to the construct at angle a and remove from self
	{
		f._i = _closest_gap(a);
		f._mo = 0; //mouse offset
		self.removeChild(f);
		fc.addChildAt(f, f._i);
		f.length = 360 * (f.meta.length / (self.len + f.meta.length));
		_set_len(self.len + f.meta.length);
		
		//_update_layout(f._i, fc.getChildAt(b(f._i - 1)).end() - (f.length / 2.0) );
		_update_layout(b(f._i-1));
	}
	
	var _leave = function(f) //remove fragment f to the construct and add to self
	{
		var g = fc.getChildIndex(f);
		f._i = undefined;
		f._mo = undefined;
		fc.removeChildAt(g);
		f.rotation = 0;
		self.addChild(f);
		_set_len(self.len - f.meta.length);
		g = b(g-1);
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
		if(n == 0) return;
		if(s==undefined) s = 0;
		if(a==undefined) a = fc.getChildAt(s).rotation;
		s = b(s);
		a = bound_angle(a);
		
		var c = a;
		var d = -1;
		var targets = new Array();
		for(var i = 0; i < n; i = i + 1)
		{
			var f = fc.getChildAt(b(i + s));
			var l = 360.0 * (f.meta.length / self.len);
			if(!f.drag)
			{
				var t = {};
				var a = false;
				if(f.rotation != c && n > 1)
				{
					a = true;
					t.rotation = c;
				}
				
				if(f.length != l)
				{
					a = true;
					t.length = l;
				}
				if(a) f.animate(t);
				t.length = l;
				t.rotation = c;
				targets.push(t);
			}
			else
			{
				d = i;
				targets.push({});
			}
			c = bound_angle(c + l);
		}
		if(d>=0 && n > 2)
		{
			var f = fc.getChildAt(b(d+s));
			var nxt = b(d+1);
			var pre = b(d-1);
			f._n = targets[nxt].rotation + targets[nxt].length / 2.0;
			f._p = targets[pre].rotation + targets[pre].length / 2.0;
		}
	}
	
	var _get_area = function(r) //return the area based on radial point r
	{
		if(r.r < radius) return CCW;
		if(r.r < radius + delta + 3 * fragw) return CW;
		return RM;
	}
	
	var _set_len = function(l) //set the length and update all the fragments
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
			s = s +'\n  ' + i + ') ' + c.meta.name + ' ['+c.meta.length+'] ';
		}
		return s;
	}
	
	var _init = function(m)
	{
		self.name = m.name;
		self.len = 0;
	//initiate the fragments
		for(var i in m.cfs)
		{
			var d = '';
			if(m.cfs[i] > 0) d = CW;
			else if (m.cfs[i] < 0) d = CCW;
			else d = undefined;
			
			fc.addChild(
				new Fragment(meta[m.cfs[i].fid],d,fragw,radius + d * delta)
			);
			self.len = self.len + meta[m.cfs[i].fid].length;
		}
		
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
		
		self.addChild(con);
		self.addChild(fc);
		
		_set_len(self.len);
		_update_layout();
		stage.update();
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
//stage is only updated if an animation is playing
var anim = 0;
var setAnim = function() {anim=anim+1;};
var clearAnim = function() {anim=anim-1; if(anim == 0) stage.update();};

var init_designer = function(cid){
	console.log('cid: ' + cid);	
	$c = $('#cdesigner');
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
		console.log('Selected: ' + m.name);
		designer.addFragment(e,m);
	}

	Ticker.setFPS(25);
	Ticker.addListener(this);
	
	stage.update();
};
var continue_init = function(m, cid)
{
	console.log('continue_init');
	library.addItems(m);
	meta = {};
	for(var i in m)
	{
		meta[m[i].fid] = m[i];
	}
	designer = new Designer(0,0,LIB_X - 2, bounds.height, cid);
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
