var canvas;
var stage;


//Declare layout values
//constants
var LIB_HANDLE_WIDTH = 20;
var LIB_STOPPER_WIDTH = 6;
var LIB_ITEM_HEIGHT = 20;
var LIB_ITEM_SPACING = 4;
//variable
var bounds = new Rectangle();
var LIB_PANEL_CLOSED = 0;
var LIB_PANEL_OPEN = 0;
var LIB_SLIDER_OPEN = 0;

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

var state = 'd';
/*
 * d - default
 * dr - drag
 * 
 * */

var library;
var slider;
var stop;
var itemPanel;
var $c;
//stage is only updated if an animation is playing
var anim = 0;
var setAnim = function() {anim=anim+1;};
var clearAnim = function() {setTimeout('anim=anim-1;', 500);};

var init_designer = function(){	
	$c = $('#cdesigner');
	$c.prop('width', $c.width());
	$c.prop('height', $c.height());
	canvas = document.getElementById('cdesigner');
	stage = new Stage(canvas);
	//stage.enableMouseOver(25);
	
	_calc_size();
	
	//Build object hirachy and set initial positions
	var hnd = _build_handle();
	stop = _build_stopper();
	stop.x = self.LIB_HANDLE_WIDTH;
	itemPanel = new Container();
	itemPanel.x = LIB_HANDLE_WIDTH + 4;
	itemPanel.y = 10;
	itemPanel.height = bounds.height - 10;
	itemPanel.width = (bounds.width - LIB_PANEL_OPEN) - LIB_HANDLE_WIDTH - LIB_STOPPER_WIDTH - 8;
	slider = new Container();
	slider.addChild(hnd);
	slider.addChild(itemPanel);
			
	library = new Container();
	library.addChild(slider);
	library.addChild(stop);
	library.x = LIB_PANEL_CLOSED;
	library.open = false;
		
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
	LIB_PANEL_CLOSED = Math.floor(bounds.width - LIB_HANDLE_WIDTH - LIB_STOPPER_WIDTH);
	LIB_PANEL_OPEN = Math.floor(0.618 * bounds.width);
	LIB_SLIDER_OPEN = LIB_PANEL_OPEN - LIB_PANEL_CLOSED;
};

var show_items = function(metadata)
{
	for(var i in metadata)
	{
		var s = new Rectangle(0, i * (LIB_ITEM_HEIGHT + LIB_ITEM_SPACING),itemPanel.width, LIB_ITEM_HEIGHT);
		var item = _draw_item(metadata[i].name, s, 'n');
		item.meta = metadata[i];
		item.onMouseOver = function(event)
		{
			console.log('item.onMouseOver()');
			_redraw_item(item, 'hl');
		}
		item.onMouseOut = function(event)
		{
			console.log('item.onMouseOut()');
			_redraw_item(item, 'n');
		}
		itemPanel.addChild(item);
	}
	stage.update();
};

//Graphics Builders
var _build_handle = function()
{
	var g = new Graphics();
	g.setStrokeStyle(1);
	g.beginFill(Graphics.getRGB(0,0,0,0.8));
	g.drawRect(LIB_HANDLE_WIDTH, 8, LIB_PANEL_OPEN, bounds.height - 16);
	
	g.beginFill(UI_BG_FILL);
	g.beginStroke(UI_BG_STROKE);
	g.drawRoundRectComplex( 0, 0, LIB_HANDLE_WIDTH, bounds.height, 4, 0, 0, 4);
	
	var s = new Shape(g);	
	
	var t = new Text('Library', '700 13.33px Lucida Grande,Lucida Sans,Arial,sans-serif', UI_FG_FILL);
	t.rotation = -90;
	t.textBaseline = 'middle';
	t.x = LIB_HANDLE_WIDTH / 2.0;
	t.y = bounds.height * 0.2; 
	
	var c = new Container();
	c.addChild(s);
	c.addChild(t);
	
	return c;
};
var _build_stopper = function()
{
	var g = new Graphics();
	g.setStrokeStyle(1);
	g.beginFill(UI_BG_FILL);
	g.beginStroke(UI_BG_STROKE);
	g.drawRoundRectComplex( 0, 0, LIB_STOPPER_WIDTH, bounds.height, 0, 4, 4, 0);
	return new Shape(g);
};
var _draw_item = function(name, size, state)
{
	var text = new Text(name, 'bold 12px arial', UI_TEXT_HL);
	var shape = new Shape(new Graphics());	
	
	text.x = size.width / 2.0;
	text.y = size.height / 2.0;
	text.textAlign = 'center';
	text.textBaseline = 'middle';
	
	var item = new Container();
	item.addChild(shape);
	item.addChild(text);
	item.state = undefined;
	item.x = size.x;
	item.y = size.y;
	item.width = size.width;
	item.height = size.height;
	return _redraw_item(item, state);
};
var _redraw_item = function(item, state)
{
	if(item.state == state) return;
	item.state = state;
	var g = item.getChildAt(0).graphics.clear();
	var text = item.getChildAt(1);
	g.setStrokeStyle(1);
	
	switch(state)
	{
		case 'hl':
			text.color = UI_TEXT_HL;
			g.beginStroke(UI_FG_STROKE_SL);
			g.beginFill(UI_FG_FILL_SL);
			_set_cursor('move');
			break;
		case 'sl':
			text.color = UI_TEXT_SL;
			g.beginStroke(UI_FG_STROKE_SL);
			g.beginFill(UI_FG_FILL_SL);
			_set_cursor('move');
			break;
		default:
			text.color = UI_TEXT;
			g.beginStroke(UI_FG_STROKE);
			g.beginFill(UI_FG_FILL);
			_set_cursor('auto');
			break;
	}
	
	g.drawRoundRect(0, 0, item.width, item.height, 4);
	return item;
}

// Mouse things --------------------------
var _on_mouse_move = function(e)
{
	if(!stage.mouseInBounds)
	{
		_closeLib();
		return;
	}
	//should the library be opened?
	if(!library.open && e.stageX > LIB_PANEL_CLOSED)
		_openLib();
	if(library.open)
	{
		if(e.stageX < LIB_PANEL_OPEN)
			_closeLib();
		else
		{
			//update states
			var item; var ht;
			for(var i = 0; i < itemPanel.getNumChildren(); i = i + 1)
			{
				item = itemPanel.getChildAt(i);
				ht = _is_over(item, e.stageX, e.stageY);
				
				if(ht && item.state != 'hl')
				{
					_redraw_item(item, 'hl');
					stage.update();
				}
				if(!ht && item.state == 'hl')
				{	
					_redraw_item(item, 'n');
					stage.update();
				}
			}
			
		}
	}
};
var _on_press = function(e)
{
	e.onMouseUp = _on_mouse_up;
	e.onMouseMove = _on_mouse_drag;
	if(library.open)
	{
		for(var i = 0; i < itemPanel.getNumChildren(); i = i + 1)
		{
			var item = itemPanel.getChildAt(i);
			var h = _is_over(item, e.stageX, e.stageY);
			if(h)
				_begin_item_drag(item);
		}
	}
}
var _on_mouse_up = function(e)
{
	i_drag = undefined;
	_set_cursor();
}
var _on_mouse_drag = function(e)
{
	if(i_drag != undefined)
	{
		i_drag.x = stage.mouseX - i_drag.width / 2.0;
		i_drag.y = stage.mouseY - i_drag.height;
		stage.update();
	}
}

var _openLib = function() 
{
	library.open = true;
	var target = {x: LIB_SLIDER_OPEN,};
	var tween = Tween.get(slider)
		.call(setAnim)
		.to(target, 500, Ease.quartOut)
		.call(clearAnim);
};
var _closeLib = function()
{
	library.open = false;
	var target = {x: 0,};
	var tween = Tween.get(slider)
		.call(setAnim)
		.to(target, 500, Ease.quartOut)
		.call(clearAnim);
};
var i_drag;
var _begin_item_drag = function(it)
{
	_closeLib();
	state = 'dr';
	i_drag = it.clone(true);
	i_drag.width = it.width;
	i_drag.height = it.height;
	i_drag.x = stage.mouseX;
	i_drag.y = stage.mouseY;
	_redraw_item(i_drag, 'sl');
	stage.addChild(i_drag);
	stage.update();
}


//Mouse Utils
var _is_over = function(ob, x,y)
{
	var l = ob.globalToLocal(x,y);
	return ob.hitTest(l.x,l.y);
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
