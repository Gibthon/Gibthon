var canvas;
var stage;


//Declare layout values
//constants
var LIB_HANDLE_WIDTH = 20;
var LIB_STOPPER_WIDTH = 6;
var LIB_ITEM_HEIGHT = 20;
var LIB_ITEM_SPACING = 4;

var FRAG_W = 25;

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
	this.dragging = false;
	this.area = 'rm'; //or ccw or cw
	
	var anim = 0; //the number of animations currently going on, controls whether to redraw
	
	this.setAnim = function(){anim = anim+1;}
	this.clearAnim = function(){anim=anim-1; if(anim==0) self.redraw();}
	
	// -------------- Transitions
	this.scale = function(s) {};
	this.rotate = function(a) {};
	this.drag = function(d) {};
	this.move = function(a) {};
	// -------------- end Transitions
	
	
	this.redraw = function()
	{
		var g = self.graphics.clear();
		switch(self.area)
		{
			case 'rm':
				g.beginFill(UI_FG_FILL)
				 .beginStroke(UI_FG_STROKE);
				 .rect(-FRAG_W * 3, -FRAG_W * .5, FRAG_W * 6, FRAG_W * .5);
				break;
			case 'cw':
				break;
			case 'ccw':
				break;
		}
	};
	this.tick = function()
	{
		if(anim > 0)
			self.redraw();
		if(self.dragging)
		{
			switch(self.area)
			{
				case 'rm':
					self.x = stage.mouseX;
					self.y = stage.mouseY;
					break;
				case 'cw':
				case 'ccw':
					
			}
		}
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

var library;

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
	library.onSelect = function(e, m)
	{
		console.log('Selected: ' + m.name);
	}
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
	if(drag == undefined) //pass to the library
		library.onMouseMove(e);
};
var _on_press = function(e)
{
	e.onMouseUp = _on_mouse_up;
	e.onMouseMove = _on_mouse_drag;
	library.onMouseDown(e);
}
var _on_mouse_up = function(e)
{
	if(drag != undefined)
	{
		_stop_dragging(drag);
		drag = undefined;
	}
	_set_cursor();
}
var _on_mouse_drag = function(e)
{
	stage.update();
}

var drag;
var _begin_item_drag = function(it)
{
	_closeLib();
	drag = new Fragment(it.meta);
	drag.drag = true;
	stage.addChild(drag);
	stage.update();
}
var _stop_dragging = function(d)
{
	switch(d.area)
	{
		case 'rm':
			stage.removeChild(d);
			d.drag = false;
			break;
		case 'cw':
			d.drag = false;
			break;
		case 'ccw':
			d.drag = false;
			break;
	}
	stage.update();
}


//Mouse Utils
var _is_over = function(ob, x,y)
{
	var l = ob.globalToLocal(x,y);
	return ob.hitTest(l.x,l.y);
}
Rectangle.prototype.containsPoint = function(p)
{
	return this.containsXY(p.x,p.y);
}
Rectangle.prototype.containsXY = function(x,y)
{
	return (x > this.x && 
			x < (this.x + this.width) &&
			y > this.y &&
			y < (this.y + this.height));
}
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

var r2s = function(r)
{
	return "("+r.x+","+r.y+")+"+r.width+"x"+r.height;
}
