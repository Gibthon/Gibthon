var canvas;
var stage;


//Declare layout values
//constants
var LIB_HANDLE_WIDTH = 20;
var LIB_STOPPER_WIDTH = 6;
//variable
var bounds = new Rectangle();
var LIB_PANEL_CLOSED = 0;
var LIB_PANEL_OPEN = 0;
var LIB_SLIDER_OPEN = 0;

var UI_BG_FILL = "rgb(100,100,100)";
var UI_BG_STROKE = "rgb(20,20,20)";
var UI_FG_FILL = "#DFEFFC";
var UI_FG_STROKE = "#C5DBEC"
var UI_TEXT = Graphics.getRGB(60,60,60);

var library;
var slider;
var stop;

var init_designer = function(){	
	var $c = $('#cdesigner');
	$c.prop('width', $c.width());
	$c.prop('height', $c.height());
	canvas = document.getElementById('cdesigner');
	stage = new Stage(canvas);
	stage.enableMouseOver();
	
	_calc_size();
	
	//Build object hirachy and set initial positions
	var hnd = _build_handle();
	stop = _build_stopper();
	stop.x = self.LIB_HANDLE_WIDTH;
	var slider = new Container();
	slider.addChild(hnd);
	//slider.addChild(itemPanel);
			
	library = new Container();
	library.addChild(slider);
	library.addChild(stop);
	library.x = LIB_PANEL_CLOSED;
	library.onMouseOver = function(event) 
	{
		var target = {x: LIB_SLIDER_OPEN,};
		var tween = Tween.get(slider)
			.to(target, 500, Ease.quartOut);
	};
	library.onMouseOut = function(event)
	{
		var target = {x: 0,};
		var tween = Tween.get(slider)
			.to(target, 500, Ease.quartOut);
	};
	
	
	stage.addChild(library);
	
	Ticker.setFPS(25);
	Ticker.addListener(this);
};
var tick = function()
{
	//console.log('tick');
	stage.update();
};
var _calc_size = function(){
	bounds.width = canvas.width;
	bounds.height = canvas.height;
	LIB_PANEL_CLOSED = bounds.width - LIB_HANDLE_WIDTH - LIB_STOPPER_WIDTH;
	LIB_PANEL_OPEN = 0.618 * bounds.width;
	LIB_SLIDER_OPEN = LIB_PANEL_OPEN - LIB_PANEL_CLOSED;
	console.log(bounds.width + 'x' + bounds.height);
};

//Graphics Builders
var _build_handle = function()
{
	var g = new Graphics();
	g.setStrokeStyle(1);
	g.beginFill(Graphics.getRGB(0,0,0,0.8));
	g.drawRect(LIB_HANDLE_WIDTH, 8, LIB_PANEL_OPEN, bounds.height - 16);
	
	g.beginFill(UI_FG_FILL);
	g.beginStroke(UI_FG_STROKE);
	g.drawRoundRectComplex( 0, 0, LIB_HANDLE_WIDTH, bounds.height, 4, 0, 0, 4);
	
	var s = new Shape(g);	
	
	var t = new Text('Library', '14px bold Arial', UI_TEXT);
	t.rotation = -90;
	t.x = LIB_HANDLE_WIDTH - 5;
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
	g.beginFill(UI_FG_FILL);
	g.beginStroke(UI_FG_STROKE);
	g.drawRoundRectComplex( 0, 0, LIB_STOPPER_WIDTH, bounds.height, 0, 4, 4, 0);
	return new Shape(g);
};
