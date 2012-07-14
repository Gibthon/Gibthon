/*
 * designer.js, Haydn King - hjking734@gmail.com - 2012
 * 
 * */

//Global things

	//the stage
	var stage;
	var $c;
	/**
	 * Animation system.
	 * */
	var _anims = 0;
	
	var setAnim = function()
	{
		_anims = _anims + 1;
	}
	var clearAnim = function()
	{
		_anims = _anims -1;
		if(_anims <= 0)
		{
			_anims = 0;
			stage.update();
		}
	}
	var tick = function()
	{
		if(_anims > 0)
			stage.update();
	}
	
	Ticker.setFPS(25);
	Ticker.addListener(this);
	
	/**
	 * Cursor
	 * */
	
	//eg auto, move, wait
	var set_cursor = function(type)
	{
		if(type == undefined) type = 'auto';
		$c.css('cursor', type);
	}

	var get_cursor = function()
	{
		return $c.css('cursor');
	}

	// prevent text cursor when dragging
	window.addEventListener("dragstart", function(fEventObject){ CancelEvent(fEventObject); } );
	window.addEventListener("selectstart", function(fEventObject){ CancelEvent(fEventObject); } );

	function CancelEvent(fEventObject) 
	{
	   if (fEventObject.preventDefault) fEventObject.preventDefault();
	   if (fEventObject.cancel != null) fEventObject.cancel = true;
	}

	/**
	 * Areas - represent the possible areas that a fragment could be in
	 * 
	**/
	var Area = {
		CCW: 0,
		CW: 1,
	};

	var AreaString = ['ccw', 'cw', 'rm',];

	var a2s = function(a) {return AreaString[a];};

//angle Functions
	var _2PI = 2 * Math.PI;
	var _D_PER_R = 360.0 / (_2PI);
	/**
	 * Convert Radians to Degrees
	 * @func r2d
	 * @param {float} r Radians
	 **/ 
	var r2d = function(r)
	{
		return r * _D_PER_R;
	};
	/**
	 * Convert Degrees to Radians
	 * @func d2r
	 * @param {float} d Degrees
	 **/ 
	var d2r = function(d)
	{
		return d / _D_PER_R;
	}
	
	/**
	 * Convert euclidian to radial (in radians)
	 * @func xy2ra
	 * @param x
	 * @param y
	 **/
	var xy2ra = function(x,y)
	{
		return {r: Math.sqrt(x*x+y*y), a: Math.atan2(y,x),};
	}
	
	/**
	 * Convert radial to euclidian
	 * @func ra2xy
	 * @param r
	 * @param a
	 **/
	var ra2xy = function(r,a)
	{
		return {x: r*Math.cos(a),y: r*Math.sin(a),};
	}

	/**
	 * Bound an angle in radians to within -PI to PI
	 * @func bound_rads
	 * @param a
	 **/
	var bound_rads = function(a) //make a within [-PI, PI]
	{
		if(a > Math.PI) return bound_rads(a - 2*Math.PI);
		if(a < -Math.PI) return bound_rads(a + 2*Math.PI);
		return a;
	};
	
	/**
	 * Bound an angle in degrees to within -180 to 180
	 * @func bound_degs
	 * @param a
	 **/
	var bound_degs = function(a) //make a within [-180, 180]
	{
		if(a > 180) return bound_degs(a - 360);
		if(a < -180) return bound_degs(a + 360);
		return a;
	};

/*
 * CanvasRenderingContext2D prototype functions
 * 
 * */
	var cp = CanvasRenderingContext2D.prototype;
	cp.moveToRA = function(r,a)
	{
		var p = ra2xy(r,a);
		return this.moveTo(p.x,p.y);
	}
	cp.lineToRA = function(r,a)
	{
		var p = ra2xy(r,a);
		return this.lineTo(p.x,p.y);
	}
	
	var g = Graphics.prototype;
	g.moveToRA = function(r,a)
	{
		var p = ra2xy(r,a);
		return this.moveTo(p.x,p.y);
	}
	g.lineToRA = function(r,a)
	{
		var p = ra2xy(r,a);
		return this.lineTo(p.x,p.y);
	}

	/**
	 * F - global fragment settings - recalculated on window change
	 * @class F
	 **/
	var F = {};
		/**
		 * Stores the radii for a each area
		 * @property radii
		 * @public
		 **/
		F.radii = [0,0,];
		/**
		 * Stores the radii for each area when a fragment is being dragged
		 * @property dradii
		 * @public
		 **/
		F.dradii = [0,0,];
		/**
		 * Stores the label radii for each area
		 * @property ldelta
		 * @public
		 **/
		F.ldelta = [0,0,];
		/**
		 * Stores the default width for a fragment
		 * @property width
		 * @public
		 **/
		F.width = 25;
		/**
		 * The length of the arrow in px
		 * @property arrow
		 * @public
		 **/
		F.arrow = 5;
		/**
		 * The maximum radius in each area
		 * @property maxRadii
		 * @public
		 **/
		F.maxRadii = [0,0,];
		F.joinRadius = 0;
		
		F.getArea = function(r)
		{
			for(var i = 0; i < 2; i = i + 1)
			{
				if(r < F.maxRadii[i])
					return i;
			}
			return undefined;
		}
		/**
		 * The minimum angle that a DisplayFragment should occupy
		 * @property minAngle
		 **/
		F.minAngle = d2r(5);
//UI
	/**
	 * The default fonts for the designer
	 * @class FONT
	 **/
	 var FONT = {
		 LG: '700 24px Lucida Grande,Lucida Sans,Arial,sans-serif',
		 M:  '500 12px Lucida Grande,Lucida Sans,Arial,sans-serif',
		 SM: '300 10px Lucida Grande,Lucida Sans,Arial,sans-serif',
	 };
	 
	 /**
	 * The default colours for the designer
	 * @class COL
	 **/
	 var COL = {
		 BLACK: 'rgb(0,0,0)',
		 WHITE: 'rgb(255,255,255)',
		 RED: 'rgb(255,0,0)',
		 GREEN: 'rgb(0,255,0)',
		 BLUE: 'rgb(0,0,255)',
	 }
	
	var _h = Math.random() * 360;
	var _grc = 0.618033988749895 * 360;
	
	var get_next_color = function()
	{
		_h = (_h + _grc) % 360;
		return Graphics.getHSL(_h,40,50);
	}

/**
* Fragment shape contains all the raw values and function to draw a fragment
* Essentially a dumb class, the FragmentShape's properties are manipulated by
* the DisplayFragment.
* @class FragmentShape
* @extends Shape
* @constructor
* @param {rotation} Rotation of the start of the fragment sits
* @param {angle} The angle corresponding to the length of the fragment
* @param {area} The initial area
* 
**/
var FragmentShape = function(_rotation, _angle, _radius)
{
	this.initialize(_rotation, _angle, _radius);
};
var fs = FragmentShape.prototype = new Shape();

//Public Properties

	/**
	 * Determines which direction the arrow should be drawn
	 * Value should be in range -1 <= clockwise <= 1
	 * Has no effect if (this.linear == true)
	 * @property clockwise
	 * @public
	 * @type float
	 **/
		fs.clockwise = 1.0;

	/**
	 * Specify the width of the fragment in px
	 * @property width
	 * @public
	 * @type int
	 **/
		fs.width = F.width;

	/**
	 * Specify the angle (length) of the fragment in rads
	 * No effect when linear == true
	 * @property angle
	 * @public
	 * @type float
	 **/
		fs.angle = 2;

	/**
	 * The radius at which to draw the fragment
	 * No effect if linear == true
	 * @property radius
	 * @public
	 * @type float
	 **/
		fs.radius = 45;

	/**
	 * The fill colour of the fragment
	 * @property fill
	 * @public
	 * @type string
	 **/
		fs.fill = Graphics.getHSL(0,0,40);
	
// constructor:
	/**
	 * @property Shape_initialize
	 * @type Function
	 * @private
	 **/
	fs.Shape_initialize = fs.initialize;
	
	/**
	 * Initialization method.
	 * @method initialize
	 * @protected
	 **/
	fs.initialize = function(_r, _a, _ra) 
	{
		this.Shape_initialize();
		this.rotation = _r;
		this.angle = _a;
		this.radius = _ra;
		this.fill = get_next_color();
	};
	
//Public Methods
	
	/**
	 * Draw a radial fragment to the context
	 * @method draw
	 * @public
	 * @param {CanvasRenderingContext2D} ctx The canvas 2D context object to draw into.
	 * @param {Boolean} ignoreCache Indicates whether the draw operation should ignore any current cache.
	 * For example, used for drawing the cache (to prevent it from simply drawing an existing cache back
	 * into itself).
	 **/
	fs.draw = function(ctx, ignoreCache)
	{
	/*	console.log('fs.draw(...)');
		p = ['width', 'angle', 'radius', 'fill'];
		for(i in p)
			console.log('	fs.'+p[i]+' = ' + this[p[i]]);
		console.log('	fs.rotation = ' + d2r(this.rotation));
	*/	
		ctx.save();
		
		ctx.lineWidth = 1;
		ctx.fillStyle = this.fill;
		ctx.strokeStyle = 'rgb(255,255,255)';
		
		ctx.beginPath();
		
		var w2 = this.width / 2.0;
		var l = this.angle;
		var a = F.arrow / this.radius;
		
		var r = this.radius;
		var d = this.clockwise;
		
		ctx.moveToRA(r - w2, 0); //left inner
		ctx.lineToRA(r, d * a); //left center
		ctx.arc(0,0,r + w2, 0, l, false); //outer edge
		ctx.lineToRA(r, l + d * a); //right center
		ctx.arc(0,0,r - w2, l, 0, true); //inner edge
				 
		ctx.fill();
		ctx.stroke();
		
		ctx.restore();
	};


/**
* 
* @class FragmentLabel
* @extends Container
* @constructor
* @param {string} The text to display
* @param {float} The radius at which to curve
* @param {float} Rotation
* 
**/
var FragmentLabel = function(_text, _radius, _rotation)
{
	this.initialize(_text, _radius, _rotation);
}
var fl = FragmentLabel.prototype = new Container();

//private variables
	/**
	 * The displayed text
	 * @property _text
	 * @private
	 * @type string
	 **/
	fl._text = '';
	
	//the angle made by the text
	fl._angle = 0.0;
	
	/**
	 * The radius at which to display - px
	 * @property _radius
	 * @private
	 * @type string
	 **/
	fl._radius = 0;
	
	/**
	 * Should the text be readable from outside the circle?
	 * @property _outward
	 * @private
	 * @type boolean
	 **/
	fl._outward = true;
	
	/**
	 * The font
	 * @property _font
	 * @private
	 * @type string
	 **/
	fl._font = FONT.SM;
	
	/**
	 * The colour
	 * @property _color
	 * @private
	 * @type string
	 **/
	fl._color = COL.BLACK;

// constructor:
	/**
	 * @property Shape_initialize
	 * @type Function
	 * @private
	 **/
	fl.Container_initialize = fl.initialize;

	//fl._cmark = null;

	/**
	 * Initialization method.
	 * @method initialize
	 * @protected
	 **/
	fl.initialize = function(_t, _ra, _ro) 
	{
		this.Container_initialize();
		
		this._text = _t;
		
		this._radius = _ra;
		
		this.rotation = r2d(_ro);
		
		//this._cmark = new Shape();
		
		this._redraw();
	};

//public methods

	/**
	 * Is the text the right way up to be read from outside the circle?
	 * @method getOutward
	 **/
	fl.getOutward = function()
	{
		return this._outward;
	};
	
	/**
	 * Is the text the right way up to be read from inside the circle?
	 * @method getInward
	 **/
	fl.getInward = function()
	{
		return !(this._outward);
	};

	/**
	 * Set whether the text can be read from outside the circle.
	 * @method setOutward
	 * @param {boolean} o outside
	 **/
	fl.setOutward = function(o)
	{
		this._outward = o;
		return this;
	};
	
	/**
	 * Set whether the text can be read from outside the circle.
	 * @method setInward
	 * @param {boolean} i inside
	 **/
	fl.setInward = function(i)
	{
		this._outward = !i;
		return this;
	};

	/**
	 * Get the font to display text
	 * @method getFont
	 **/
	fl.getFont = function()
	{
		return this._font;
	};
	
	/**
	 * Set the font to display text
	 * @method setFont
	 * @param {string} _font The new font
	 **/
	fl.setFont = function(_font)
	{
		this._font = font;
		return this;
	};
	
	/**
	 * Get the colour to display text
	 * @method getColor
	 **/
	fl.getColor = function()
	{
		return this._color;
	};
	
	/**
	 * Set the color to display text
	 * @method setColor
	 * @param {string} _color The new colour
	 **/
	fl.setColor = function(_color)
	{
		this._color = color;
		return this;
	};
	
	/**
	 * Get the displayed text
	 * @method getText
	 **/
	fl.getText = function()
	{
		return this._text;
	};
	
	/**
	 * Set the displayed text
	 * @method setText
	 * @param {string} _text The text to set
	 * @param {boolean} _update Whether to update the stage
	 * @default true
	 **/
	fl.setText = function(_text, _update)
	{	
		this.text = _text;
		this._redraw();
		//Cache
	};
	
	/**
	 * Get the radius
	 * @method getRadius
	 **/
	fl.getRadius = function()
	{
		return this._radius;
	};
	
	/**
	 * Set the radius
	 * @method setRadius
	 * @param {float} _radius The new radius
	 * @param {boolean} _update Whether to update the stage
	 * @default true
	 **/
	fl.setRadius = function(_radius)
	{		
		this._radius = _radius;
		this._realign();
	}
	
	/**
	 * Get the rotation
	 * @method getRadius
	 **/
	fl.getRotation = function()
	{
		return deg2rad(this.rotation) + 0.5 * this._angle;
	};
	
	/**
	 * Set the rotation
	 * @method setRotation
	 * @param {float} _rotation The new rotation (rads)
	 * @param {boolean} _update Whether to update the stage
	 * @default true
	 **/
	fl.setRotation = function(_rotation)
	{	
		//find overall rotation
		var r = bound_degs(r2d(_rotation) + this.parent.rotation);
		//nb. 	r < 0 => outward = false
		//		r > 0 => outward = true 
		//if the outwardness needs changing
		if( (r > 0) != this._outward)
		{
			this._outward = !this._outward;
			this._realign();
		}
		
		this.rotation = r2d(_rotation) - r2d(0.5 * this._angle);
		return this;
	};
	
	/**
	 * Get the visibility
	 * @method isVisible
	 **/
	fl.isVisible = function()
	{
		return this.visible;
	};
	
	/**
	 * Show the label
	 * @method show
	 **/
	fl.show = function()
	{
		this.visible = true;
		return this;
	};
	
	/**
	 * Hide the label
	 * @method hide
	 **/
	fl.hide = function()
	{
		this.visible = false;
		return this;
	};
	
	fl.toString = function()
	{
		return '[ FragmentLabel (text='+this._text+')]';
	}
	
//private methods
	/**
	 * Redraw all the letters
	 * @method _redraw
	 * @protected
	 **/
	fl._redraw = function()
	{
		//clear any previous letters
		this.removeAllChildren();
	
		//add the new letters, one at a time
		for(var i in this._text)
		{
			var l = new Text(this._text.charAt(i), this._font, this._color);
			this.addChild(l);
		}
		
		//this.addChild(this._cmark);
		return this._realign();
	};
	
	/**
	 * Redraw all the letters
	 * @method _redraw
	 * @protected
	 **/
	fl._realign = function()
	{		
		var r = 0;
		var erad = this._radius;
		if(this._outward)
			erad = - this._radius;
		
		//for each letter
		for(var i = 0; i < this._text.length; i = i + 1)
		{
			var l = this.getChildAt(i);
			//reset the center
			l.x = 0;
			l.maxWidth = 1000;
			//set the y-offset accordingly
			l.regY = erad;
				
			//set the rotation
			l.rotation = r2d(r) + 90;
			if(this._outward)
			{
				l.textBaseline = 'top';
				l.rotation = -r2d(r) - 90;
			}
			else
				l.textBaseline = 'bottom';
			
			r = r + l.getMeasuredWidth() / this._radius;
			
			
		}
		//set the rotation
		if(this._outward) r = - Math.abs(r);
		else r = Math.abs(r);
		
    //console.log('this.rotation = r2d((d2r('+this.rotation+') + 0.5 * '+this._angle+') - 0.5 * '+r+');');
		this.rotation = r2d((d2r(this.rotation) + 0.5 * this._angle) - 0.5 * r);
    //console.log('  = '+this.rotation);
		this._angle = r;
		
		//cache the result
		var pad = 20;
		var dr = 0;
		//this._outward ? dr = 90 : dr = -90;
		p1 = ra2xy(this._radius, 0);
		p2 = ra2xy(this._radius, this._angle);
		
		var x = Math.min(p1.x, p2.x) - pad;
		var y = Math.min(p1.y, p2.y) - pad;
		var w = Math.abs(p1.x - p2.x) + 2*pad;
		var h = Math.abs(p1.y - p2.y) + 2*pad;
		
		
	/*	cache debugging
	 * console.log(this+'.cache('+x+','+y+','+w+','+h+');');
	 * this._cmark.graphics
			.clear()
			.setStrokeStyle(1)
			.beginStroke(COL.RED)
			.rect(x,y,w,h)
			.endStroke()
			.closePath()
			.beginFill(COL.RED)
			.drawCircle(p2.x,p2.y,3)
			.beginFill(COL.GREEN)
			.drawCircle(p1.x,p1.y,3);
	*/	
		
		this.cache(x,y,w,h);
		
		return this;
	};
	
	fl.toString = function() {return '[ FragmentLabel ("'+this._text+'")]';};
	
/**
* 
* @class DisplayFragment
* @extends Container
* @constructor
* @param {Fragment} _f The Fragment which this represents
* @param {ConstructFragment} _cf The ConstructFragment - represents it's position in this construct
* @default null
* 
**/
var DisplayFragment = function(_f, _cf)
{
	this.initialize(_f, _cf);
}
var df = DisplayFragment.prototype = new Container();

//Public Variables

//Private Variables
	 
	/**
	 * The Fragment Object
	 * @property _f
	 * @private
	 * @type Fragment
	 **/
	 df._f = null;
	 
	 /**
	 * The ConstructFragment
	 * @property _cf
	 * @private
	 * @type ConstuctFragmnet
	 **/
	 df._cf = null;
	 
	 /**
	 * The FragmentShape draws the fragment to the canvas
	 * @property _fs
	 * @private
	 * @type FragmentShape
	 **/
	 df._fs = null;
	 
	 /**
	  * store a copy the the fragment shape's properties to avoid showing tweend values to the outside world
	  * */
	 df._props = null;
	 
	 /**
	 * The FragmentLabel draws the label to the canvas
	 * @property _fl
	 * @private
	 * @type FragmentLabel
	 **/
	 df._fl = null;
	 
	 /**
	 * Are we being dragged?
	 * @property _drag
	 * @private
	 * @type boolean
	 **/
	 df._drag = null;
	 
	 /**
	 * Which area are we in?
	 * @property _area
	 * @private
	 * @type Area
	 **/
	 df._area = null;
	
	/**
	 * The angle between the start and the postion at which the mouse was clicked
	 * @property _mouse_offset
	 * @private
	 * @type float
	 **/
	 df._mouse_offset = 0.0;
	 
	 df._mouse_down = new Point(0,0);
	 
	 df._mousedownEvent = null;
	 
//Constructor
	/**
	 * @property Shape_initialize
	 * @type Function
	 * @private
	 **/
	df.Container_initialize = df.initialize;

	/**
	 * Initialization method.
	 * @method initialize
	 * @protected
	 **/
	df.initialize = function(_f, _cf) 
	{
		this.Container_initialize();
		//set initial values
		
		this._props = {
			'rotation': 0.0,
			'angle':0.0,
			'radius':0.0,
			'alpha':1.0,
			'clockwise':1.0,	
		};
		
		this._f = _f;
		this._cf = _cf;
		this._area = Area.CW;
		if(_cf != undefined)
		{
			if(_cf.strand == 1)
				this._area = Area.CW;
			if(_cf.strand == -1)
				this._area = Area.CCW;
		}

		this._fl = new FragmentLabel(
            this._f.getName(), 
            F.radii[this._area] + F.ldelta[this._area], 
            this._props.angle / 2.0);

		this._fs = new FragmentShape(0,0,F.radii[this._area]);
		this._props.radius = F.radii[this._area];
		
		this._drag = false;
		
		if(this._area == Area.CCW)
		{
			this.setClockwise(-1);
		}
				
		this.addChild(this._fs);
		this.addChild(this._fl);
		
	};
//Public Methods
//Getters & Setters
	/**
	 * Get the current area
	 * @method getArea
	 **/
	df.getArea = function()
	{
		return this._area;
	}
	 
	 /**
	 * Set the current area
	 * @method setArea
	 * @param {Area} _a The new area
	 **/
	df.setArea = function(_a)
	{
		this.x = 0; this.y = 0;
		this._area = _a;
        if(this._cf)
            this._cf.strand = (this._area == Area.CW) ? 1 : -1;
		this._drag ?
			this.setRadius(F.dradii[_a]) :
			this.setRadius(F.radii[_a]);
		
		(_a == Area.CW) ?
			this.setClockwise(1.0) :
			this.setClockwise(-1.0);
		
		this._setLabelRadius();
	}
	
	df.setRotation = function(rads)
	{
		this._props.rotation = this._fs.rotation = r2d(bound_rads(rads));
        this._fl.setRotation(this.getMid());
	}
	
	df.setAngle = function(angle)
	{
		this._props.angle = this._fs.angle = angle;
        this._fl.setRotation(this.getMid());
	}
	
	df.setRadius = function(px)
	{
		this._props.radius = this._fs.radius = px;
	}
	
	df.setAlpha = function(a)
	{
		this._props.alpha = this._fs.alpha = a;
	}
	
	df.setClockwise = function(cw)
	{
		this._props.clockwise = this._fs.clockwise = cw;
	}
	
	df.getRotation = function(rads)
	{
		return d2r(this._props.rotation);
	}
	
	df.getMouseRotation = function(rads)
	{
		return d2r(this._props.rotation) + this._mouse_offset;
	}
	
	df.getAngle = function(angle)
	{
		return this._props.angle;
	}
	
	df.getRadius = function(px)
	{
		return this._props.radius;
	}
	
	df.getAlpha = function(a)
	{
		return this._props.alpha;
	}
	
	df.getClockwise = function(cw)
	{
		return this._props.clockwise;
	}
	
	df.f = function() {return this._f;};
	
	/**
	 * Return whether the fragment is being dragged or not
	 * @method isDragging
	 **/
	df.isDragging = function()
	{
		return this._drag;
	}
	
	/**
	 * Return the current total length of the fragment (bp)
	 * @method getLength
	 **/
	df.getLength = function()
	{
		if(this._cf != undefined)
			return this._cf.getLength();
		return this._f.getLength();
	}
	
	/**
	 * Return the current order of the fragment, -1 if not in the construct
	 * @method getOrder
	 **/
	df.getOrder = function()
	{
		if(this._cf != null)
			return this._cf.order;
		return -1;
	}
	
	df.getCf = function()
	{
		return this._cf;
	}

    df.getServer = function()
    {
        return this.parent.getServer();
    }
	
	/**
	 * Get the angle (rads) at which the fragment starts
	 * @method getStart
	 **/
	df.getStart = function()
	{
		return d2r(this._props.rotation);
	}
	
	/**
	 * Get the angle (rads) at the middle
	 * @method getMid
	 **/
	df.getMid = function()
	{
		return d2r(this._props.rotation) + this._props.angle / 2.0;
	}
	
	/**
	 * Get the angle (rads) at which the fragment ends
	 * @method getEnd
	 **/
	df.getEnd = function()
	{
		return d2r(this._props.rotation) + this._props.angle;
	}
	
	/**
	 * Get the angle (rads) of the fragment (i.e. length)
	 * @method getAngle
	 **/
	df.getAngle = function()
	{
		return this._props.angle;
	}
	
	/**
	 * Animate the properties of the framgment
	 * @method animate
	 * @param {Object} t The target properties
	 * Supports rotation (rads), angle(rads), radius, alpha, clockwise
	 * @param {function} cb A callback for when the animation completes 
	 **/
	df.animate = function(t, cb)
	{
        var change = false;
		var hide = false; //should the label be hidden?
		var self = this;
		
	
	//check for changes
		var delta = 1 / this._props.radius;
		//make sure any rotation goes the fastest direction
		if(t.rotation != undefined)
		{
			var r = bound_degs(r2d(t.rotation) - this._props.rotation); //r is the distance and direction of shortest rotation
			t.rotation = this._props.rotation + r;
			if(Math.abs(r) > delta) change = true;
            hide = true;
		}
		//angle
		if(t.angle != undefined)
		{
			if(Math.abs(this._props.angle - t.angle) > delta) change = true;
		}
		//radius
		if(t.radius != undefined)
		{
			if(Math.abs(this._props.radius - t.radius) > 1.0) 
			{
				change = true;
				hide = true;
			}
		}
		//alpha
		if(t.alpha != undefined)
		{
			if(this._props.alpha != t.alpha) change = true;
		}
		//clockwise
		if(t.clockwise != undefined)
		{
			if(this._props.clockwise != t.clockwise) change = true;
		}
		
		if(!change) return;

		//update the saved values
		this._set(t);
		this._props.rotation = bound_degs(this._props.rotation);
		
		var done = function() {
            self.setRotation(d2r(self._fs.rotation));
            self._setLabelRadius();
            self._fl.show();
		};
		if(hide) this._fl.hide();
		
		var tween = Tween.get(this._fs)
		 .call(setAnim) //enable global animation
		 .to(t, 250, Ease.quartOut)
		 .call(done)
		 .call(clearAnim); //disable global animation
		if(cb != undefined) tween.call(cb);
        stage.update();
	}
	
	/**
	 * Set the properties of the framgment
	 * @method set
	 * @param {Object} t The target properties
	 * Supports start (rads), angle(rads), radius, alpha, clockwise
	 **/
	df.set = function(t)
	{
		//set rotation
		if(t.rotation != undefined)
		{
			this.setRotation(t.rotation);
		}
		//set angle
		if(t.angle != undefined)
		{
			this.setAngle(t.angle);
		}
		//set radius
		if(t.radius != undefined)
		{
			this.setRadius(t.radius);
		}
		//set alpha
		if(t.alpha != undefined)
		{
			this.setAlpha(t.alpha);
		}
		//set clockwise
		if(t.clockwise != undefined)
		{
			this.setClockwise(t.clockwise);
		}
	}
	
	df.toString = function()
	{
		return '[DisplayFragment (f.name = "'+this._f.getName()+'")]'; 
	}
	
//Mouse Events
	/**
	 * Handle a mouse Over
	 * @method onMouseOver
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onMouseOver = function(ev)
	{
		if(!this._drag && (get_cursor() == 'auto'))
		{
			set_cursor('pointer');
			this.alpha = 0.8;
			stage.update();
		}
	}
	
	/**
	 * Handle a mouse Out
	 * @method onMouseOut
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onMouseOut = function(ev)
	{
		if(!this._drag && (get_cursor() == 'pointer'))
		{
			set_cursor();
			this.alpha = 1.0;
			stage.update();
		}
	}
	
	/**
	 * Handle a mouse click
	 * @method onClick
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onClick = function(ev)
	{
		//if we were about to start dragging
		if(!this._drag)
		{
			//do on-click stuff
			this.parent.designer().showInfo(this);
			return;
		}
		this._drag = false;
		this._mousedownEvent = null;
	}
	
	/**
	 * Handle a mouse press
	 * @method onPress
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onPress = function(ev)
	{
		var self = this;
		
		this._mousedownEvent = ev;
		this._mouse_down = this._get_mev(ev);
		
		this._mouse_offset = bound_rads(this._mouse_down.a - d2r(this._props.rotation));
		stage.onMouseMove = function(e) {self.onDrag(e);};
		stage.onMouseUp = function(e) {stage.onMouseMove = null; stage.onMouseUp = null;};
	}
	
	df.onDragStart = function(ev)
	{
		this._drag = true;
		this.setRadius(F.dradii[this._area]);
		this._setLabelRadius();
		
		set_cursor('move');
		
		stage.update();
		
		var self = this;
		//listen for mouse ups
		stage.onMouseMove = function(e) {self.onDrag(e);};
		stage.onMouseUp = function(e) {
            self.onDrop(e);
        };
	}
	
	/**
	 * Handle a mouse Drag
	 * @method onDrag
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onDrag = function(ev)
	{
		var p = this._get_mev(ev);
		if(this._drag)
		{
			this.setRotation(p.a - this._mouse_offset);
			
			var area = F.getArea(p.r);
			if(area != this._area)
			{
				if(area != undefined)
				{
					this.setArea(area);
				}
				else
				{
					stage.onMouseMove = null;
					stage.onMouseUp = null;
					ev.onMouseMove = null;
					
					set_cursor();
					this.parent.designer().leave(this);
					return;
				}
			}
			
			
			stage.update();
			this.parent.sortOne(this);
		}
		else
		{
			var x = (p.x - this._mouse_down.x);
			var y = (p.y - this._mouse_down.y);
			if((x*x + y*y) > (0.5 * F.width * F.width))
				this.onDragStart(ev);
		}
		
	};
	
	/**
	 * Handle a mouse drop
	 * @method onDrop
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onDrop = function(ev)
	{
        var self = this;
		stage.onMouseMove = null;
		stage.onMouseUp = null;
		
        t = {radius: F.radii[this._area],};
		t.rotation = this.parent.onDrop(this);
		this._mousedownEvent = null;
	    set_cursor();
		this.animate(t, function() 
        {
            //stop dragging when the animation finishes to prevent erroneous
            //clicks from being generated
            self._drag = false;
        });

        if(this._cf == undefined)
            {
                this.getServer().addFrag(this);
            }
        else
            {
                this.getServer().reorder(this.parent.children);
            }
	}
	

//Private Methods
	
	df._setLabelRadius = function()
	{
		var r = this.getRadius() + F.ldelta[this._area];
		this._fl.setRadius(r);
	}
	
	/**
	 * Get the local radial and euclidian mouse position
	 * @method _get_mev
	 * @private
	 * @param {MouseEvent} ev The MouseEvent
	 **/
	df._get_mev = function(ev)
	{
		var p = this.globalToLocal(ev.stageX, ev.stageY);
		var rp = xy2ra(p.x,p.y);
		return { x: p.x, y:p.y, r: rp.r, a: rp.a,};
	}
	
	/**
	 * Set the prop values to t
	 * */
	df._set = function(t)
	{
        for(var i in this._props)
		{
			if(t[i]!=undefined)
			{
				this._props[i] = t[i];
			}
		}
	}


/**
* Manitains the position of the fragments while they're in the construct
* @class FragmentContainer
* @extends Container
* @constructor
* @param {Server} server A server object for saving changes
**/
var FragmentContainer = function(designer)
{
	this.initialize(designer);
}
var fc = FragmentContainer.prototype = new Container();

//private Variables
	
	/**
	 * @property _length
	 * @private
	 **/
	fc._length = 0;
	
	/**
	 * Effective length - used to calculate fragment angles
	 * @property _eff_length
	 * @private
	 **/
	fc._eff_length = 0;
	
	/**
	 * The minimum length that a display fragment should be set
	 * @property _lmin
	 * @private
	 **/
	fc._lmin = 0;
	
	/**
	 * Length to be added for a gap
	 * @property _gap_length
	 * @private
	 **/
	fc._gap_length = 0;
	
	fc._designer = null;

//Constructor
	/**
	 * @property Container_initialize
	 * @type Function
	 * @private
	 **/
	fc.Container_initialize = fc.initialize;

	/**
	 * Initialization method.
	 * @method initialize
	 * @protected
	 **/
	fc.initialize = function(d) 
	{
		this.Container_initialize();
		this._designer = d;	
	};
	
//public methods

	fc.designer = function() {return this._designer;};
	/**
	 * Add a new DisplayFragment
	 * @method addFragAt
	 * @param {DisplayFragment} df The DisplayFragment to add
	 * @param {int} pos The position to add the DisplayFragment
	 * @default df.cf.order
	 **/
	fc.addFragAt = function(df, pos)
	{
		if(!$.isNumeric(pos))
			(df._cf == undefined) ? pos = 0 : pos = df._cf.order;
		else
			pos = this._bound(pos);
		
		this.addChildAt(df, pos);
		
		this._updateLength();
		
		df._mouse_offset = Math.PI*df.getLength()/this._eff_length;
		this._updateLayout( 
			pos+1,
			this.getFragAt(pos+1).getStart() + df._mouse_offset,
			true
		);
		
		return this;
	}

    fc.getServer = function()
    {
        return this._designer.getServer();
    }
	
    fc.debug = function()
    {
        console.log('----------------------------------------------\n'+
                    'FragmentContainer debug information:');
        //console.log('  ._length = '+this._length);
        console.log('  ._eff_length = '+this._eff_length);
        //console.log('  numChildren: '+this.getNumChildren());
        for(var i = 0; i<this.getNumChildren(); i=i+1)
        {
            var c = this.getFragAt(i);
            console.log('  ['+i+'] - '+c+' '+c._f.getLength()+'bp '+
                        c.getStart()+' -> '+c.getEnd()+' _drag = '+c._drag);
            //console.log('  ._props.rotation: '+c._props.rotation);
            //console.log('  ._fc.rotation: '+c._fs.rotation);
        }
        console.log('----------------------------------------------');
    }

	fc.add = function(df)
	{
		//figure where the fragment should be added
		//mouse position
		var xy = this.globalToLocal(stage.mouseX, stage.mouseY);
		var p = xy2ra(xy.x,xy.y);


		var min_v = _2PI; var min_i = 0;
		for(var i = 0; i < this.children.length; i = i + 1)
		{
			var v = Math.abs(bound_rads(p.a - this.children[i].getStart()));
			if(v < min_v)
			{
				min_v = v;
				min_i = i;
			}
		}
		//add the fragment
		
		df._drag = true;
		this.addFragAt(df, min_i);
		df.setRotation(p.a - df._mouse_offset);
		df.setAngle(_2PI*df.getLength()/this._eff_length);

		//begin the dragging!
		df.onDragStart($.Event('mousemove'));
	}
	
	/**
	 * Add multiple fragments (e.g. at initial load)
	 * @method addMulti
	 * @param {Array[DisplayFragment]} dfs
	 **/
	fc.addMulti = function(dfs)
	{
		for(i in dfs)
			this.addChild(dfs[i]);
		
        this.children.sort(function(a,b){
            return a.getCf().order - b.getCf().order;
        });
		
		this._updateLength();
		this._updateLayout(0,0,false);
		
		return this;
	}
	
	/**
	 * Remove a displayFragment
	 * @method rm
	 * @param {DisplayFragment} df The DisplayFragment to remove
	 **/
	fc.rm = function(df)
	{
		var i = this.getChildIndex(df);
		if(i < 0) return this;
		
		var a = df.getMid();
		
		this.removeChildAt(i);
		this._updateLength();
		this._updateLayout(this._bound(i), a);
		this._datum = NaN;
		return this;
	}
	
	/**
	 * Notify the container that a drop occured
	 * @method drop
	 * @param {DisplayFragment} df
	 **/
	fc.onDrop = function(df)
	{
		this._datum = NaN;
		return this.getFragAt(this.getChildIndex(df)-1).getEnd();
	}
	
	/**
	 * Get the length of the construct in base pairs
	 * @method getLength
	 **/
	fc.getLength = function()
	{
		return this._length;
	}
	
	/**
	 * get a DisplayFragment by index, with wrapping
	 * @method getFragAt
	 * @param {int} i The Offset
	 **/
	fc.getFragAt = function(i)
	{
		return this.getChildAt( this._bound(i));
	}
	
	/**
	 * Get a DisplayFragment by cfid, or undefined if it doesn't exist
	 * @method getFragByCFID
	 * @param {int} cfid The ConstructFragment ID to search for
	 **/
	fc.getFragByCFID = function(cfid)
	{
		for(var i = 0; i < this.getNumChildren(); i = i + 1)
		{
			var df = this.getChildAt(i);
			if(df.cf != undefined)
			{
				if(df.cf.id == cfid)
					return df;
			}
		}
		return undefined;
	}
	
	
//PRIVATE METHODS
	/**
	 * Update the layout
	 * @method _updateLayout
	 * @param {int} startFrag
	 * @default 0
	 * @param {float} startAngle
	 * @default startFrag.rotation
	 * @param {bool} animate Whether the layout change should be animated
	 * @default true
	 * @private
	 **/
	fc._updateLayout = function(startFrag, startAngle, animate)
	{
		//sort out the values
		if(startFrag == undefined)
			startFrag = 0;
		else
			startFrag = this._bound(parseInt(startFrag));
			
		if(startAngle == undefined)
			startAngle = this.getChildAt(startFrag).getStart();
		else
			startAngle = bound_rads(startAngle);
		
		if(animate == undefined)
			animate = true;
		
		//keep track of the current rotation
		var r = startAngle;
		//keep hold of the targets
		var targets = new Array();
		
		//for each fragment
		for(var i = 0; i < this.getNumChildren(); i = i + 1)
		{
			//get the appropriate fragment
			var df = this.getFragAt(i + startFrag);
			
			//make the target
			var t = {};
			
			if(!df.isDragging())
			{
				t.rotation = bound_rads(r);
			}
			
			//calculate the angle
			var l = df.getLength();
			if(l < this._lmin) l = this._lmin;
			
			var a = _2PI *  l / this._eff_length;
			
			t.angle = a;
			r = r + a;
			targets.push(t);
		}

		//apply the targets all at once
		if(animate)
		{
			for(var t = 0; t < targets.length; t = t+1)
			{
				this.getFragAt(t + startFrag).animate(targets[t]);
			}
		}
		else
		{
			for(var t = 0; t < targets.length; t = t+1)
			{
				this.getFragAt(t + startFrag).set(targets[t]);
			}
		}
		
		return this;
	}
	
	/**
	 * sort the fragment s, assuming all others are in order
	 * @method sortOne
	 * @param {DisplayFragment} s
	 * @public
	 **/
	
	fc._datum = NaN; 
	fc._lim = NaN;
	
	//fc._gdb = new Shape();
	
	fc.sortOne = function(s)
	{
		if(isNaN(this._datum))
		{
			var i = this.getChildIndex(s);
			var n = this.getFragAt(i+1);
			var p = this.getFragAt(i-1);
			
			this._lim = 0.5 * (0.5*n.getAngle() + 0.5*p.getAngle() + s.getAngle());
			this._datum = bound_rads(n.getMid() - this._lim);
			
			if(!this.parent.contains(this._gdb))
				this.parent.addChild(this._gdb);
			
			/*this._gdb.graphics.clear()
			.moveTo(0,0)
			.beginStroke(COL.GREEN)
			.lineToRA(100, this._datum)
			.endStroke()
			.moveTo(0,0)
			.beginStroke(COL.RED)
			.lineToRA(120, this._datum + this._lim)
			.endStroke()
			.moveTo(0,0)
			.beginStroke(COL.RED)
			.lineToRA(120, this._datum - this._lim)
			.endStroke();*/
			
		}
		
		//test for compliance
		var a = bound_rads(s.getMouseRotation() - this._datum);
		if(a > this._lim)
		{
			this._swap(this.getChildIndex(s), 1);
			this._datum = NaN;
		}
		if(a < -this._lim)
		{
			this._swap(this.getChildIndex(s), -1);
			this._datum = NaN;
		}
		
		return this;
	}
	
	//a = offset to swap, offset = 1 | -1
	fc._swap = function(i,offset)
	{
		var ipo = this._bound(i + offset);
		var imo = this._bound(i - offset);
		
		var t = this.children[i];
		this.children[i] = this.children[ipo];
		this.children[ipo] = t;
		
		var target = {};
		switch(offset)
		{
			case -1:
				target.rotation = this.children[imo].getStart() - this.children[i].getAngle();
				break;
			case 1:
				target.rotation = this.children[imo].getEnd();
				break;
		}
		this.children[i].animate(target);
	}
	
	
	/**
	 * sort the fragments by rotation
	 * @method _sortAll
	 * @public
	 **/
	fc.sortAll = function()
	{
		var s = function(a,b)
		{
			return bound_rads(a.getStart() - b.getStart());
		}
		children.sort(s);
		_updateLayout();
		return this;
	}
	
	/**
	 * Update the calculated length, eff_length and gap_length, tell the designer if it changed
	 * @method _updateLength
	 * @private
	 **/
	fc._updateLength = function()
	{
		var _old = this._length;
		
		//find the new real length, and store all the lengths
		this._length = 0;
		var lengths = new Array();
		for(var i = 0; i < this.getNumChildren(); i = i + 1)
		{
			var l = this.getChildAt(i).getLength();
			this._length = this._length + l;
			lengths.push(l);
		}
		
        this._designer.setLength(this._length);
        
		//find the minimum length
		this._lmin = Math.ceil(this._length * F.minAngle / (_2PI));
		
		//find the effective length
		this._eff_length = 0;
		for(var i in lengths)
		{
			if(lengths[i] > this._lmin)
				this._eff_length = this._eff_length + lengths[i];
			else
				this._eff_length = this._eff_length + this._lmin;
		}
		return this;
	}
	
	/**
	 * Return i bounded to the number of children
	 * @method _bound
	 * @param {int} i
	 **/
	fc._bound = function(i)
	{
        if(this.children.length == 0) return 0;
		i = i % this.children.length;
		if(i < 0) return i + this.children.length;
		return i;
	}
	
	fc.toString = function() {return '[FragmentContainer (numChildren='+this.getNumChildren()+')]';};

	
/**
* Handles server activity and displays status messages
* @class Server
* @extends Container
**/
var Server = function()
{
	this.initialize();
}
var s = Server.prototype = new Container();

	//private data
	
	s._text = new Text('', FONT.M, COL.BLACK);
	s._sm_text = new Text('', FONT.SM, COL.RED);
    
    s._con = null;


	s._container_init = s.initialize;

	s.initialize = function()
	{
		this._text.x = 0;
		this._text.maxWidth = 1000;
		this._sm_text.x = 10;
		this._sm_text.y = -this._text.getMeasuredLineHeight() - 2;
		this._sm_text.visible = false;
		this._sm_text.maxWidth = 1000;
		this._text.visible = false;
		
		this.addChild(this._text, this._sm_text);
	}

	//Public Functions

    s.getConstruct = function()
    {
        return this._con;
    }
	s.getInfo = function(cid, cb)
	{
		this._showMessage('Getting Info...');
		
		var self = this;
        libDesigner.getConstructByID(cid, function(c)
        {
             self._con = c;
             self._hideMessage();
             if($.isFunction(cb)) cb(c);
        });
		return this;
	}
	
	s.addFrag = function(df)
	{
		this._showMessage('Adding Fragment...');
        var d = (df._area == Area.CW) ? 1 : -1;
		var self = this;
		this._con.addFragment(df._f, df.parent.getChildIndex(this), d, 
            function(cf)
			{
				self._hideMessage();
                console.log('setting _cf = '+cf);
                df._cf = cf;
			});
	}
	
	s.rmFrag = function(cfid, cb)
	{
		this._showMessage('Removing Fragment...');
		var self = this;
		this._con.rmFragment(cfid, function(data)
			{
				self._hideMessage();
				if($.isFunction(cb)) cb();
			});
	}

    s.reorder = function(dfs, cb)
    {
        this._showMessage('Reordering Fragments...');
        var self = this;
        var cfids = [];
        var dirs = [];
        for(var i in dfs)
            {
                cfids.push(dfs[i]._cf.id);
                console.log('dfs['+i+']._cf.strand: '+dfs[i]._cf.strand);
                dirs.push(dfs[i]._cf.strand);
            }
        this._con.reorder(cfids,dirs,function()
          {
              self._hideMessage();
              if($.isFunction(cb)) cb();
          });
    }

	//Private functions
	
	s._handle_error = function(desc, textStatus, errorThrown)
	{
		this._showMessage(	"Sorry, there was an error while " + desc + ", please try reloading.",
							"status: '" + textStatus + "' error: '" + errorThrown + '"');
	}
	
	s._handle_success = function(cb, data)
	{
		this.visible = false;
		this._sm_text.visible = false;
		cb(data);
	}
	
	s._showMessage = function(msg, emsg)
	{
		this._text.text = msg;
		this._text.color = COL.BLACK;
		this._text.visible = true;
		if(emsg)
		{
			s._sm_text._text = emsg;
			s._sm_text.visible = true;
			this._text.color = COL.RED;
			this._sm_text.color = COL.BLACK;
		}
		stage.update();
	}
	
	s._hideMessage = function()
	{
		s._text.visible = false;
		s._sm_text.visible = false;
		stage.update();
	}
	
	s.toString = function() {return '[Server ]'};

/**
* The actual designer, called with a jQuery HTML canvas
* @class Designer
* @extends Container
**/
var Designer = function($canvas, cid)
{
	this.initialize($canvas, cid);
}
var d = Designer.prototype = new Container();
	
	d._$canvas = null;
	d._$info = null;
	
	d._cid = 0;
	
	d._width = 0;
	d._height = 0;
	
	//containers
	d._child = null;
	d._server = null;
	d._fc = null;
	d._tname = null;
	d._tlen = null;

	//Constuctor
	d._container_initialize = d.initialize;
	d.initialize = function($canvas, cid)
	{
		var self = this;
		this._$canvas = $canvas;
		$c = $canvas;
		this._cid = cid;
		$(window).resize(function() {self._calcSize(); stage.update();});
		
		this._child = new Container();
		this._tname = new Text('', FONT.LG, COL.BLACK);
		this._tlen = new Text('', FONT.M, COL.BLACK);
		this._tname.textAlign = 'center';
		this._tname.maxWidth = 1000;
		this._tlen.textAlign = 'center';
		this._tlen.maxWidth = 1000;
		
		this._fc = new FragmentContainer(this);
//$(window).keypress(function() {self._fc.debug();});
		
		this._server = new Server();
		
		this._child.addChild(this._tname, this._tlen, this._fc);
		this.addChild(this._child, this._server);
		
		this._calcSize();
		
		this.setName('My Name (' + this._cid + ')');
		this.setLength(150);
		
		//setup the canvas
		stage = new Stage(document.getElementById('cdesigner')); //$canvas.get());
		stage.addChild(this);
		stage.enableMouseOver(15);
		
		var self = this;
		self._server.getInfo(self._cid, function(c) {self._gotInfo(c)});
		
		self._initInfo();
		
		stage.update();

        //listen for items being dragged into the canvas
		var o = this._$canvas.parent().offset();
        $canvas.droppable({
            accept: '.jFragment',
            over: function(event, ui) {
                var jf = ui.draggable;
                jf.on('drag', function(event, ui) {
                    var p = self._fc.globalToLocal( event.pageX - o.left, 
                                                   event.pageY - o.top);
                    if( (p.x*p.x + p.y*p.y) < F.joinRadius * F.joinRadius )
                    {
                        //try and avoid being clobbered by high speed
                        //mice
                        jf.off('drag');
                        if(jf.data('cf'))
                        {
                            jf.off('dragstop');
                            jf.on('dragstop', function(ev,ui){
                                jf.remove(); //remove from the DOM
                                //but don't tell the server anything
                            });
                        }
                        self.join(jf);
                        //stop ui.mouse from getting confused...
                        $(document).mouseup();
                        return false; //triggers 'stop'
                    }
                    return true;
                });
            },
        });


        //debug
        $canvas.on('mouseup', function(){ console.log('canvas mouseup');});

	}
	
    d.getConstruct = function()
    {
        return this._server.getConstruct();
    }
    d.updateStage = function()
    {
        stage.update();
    }
	d.getName = function()
	{
		return this._tname.text;
	}
	d.setName = function(name)
	{
		this._tname.text = name;
	}
	d.setLength = function(len)
	{
		this._tlen.text = len + ' bp';
	}
	d.getLength = function()
	{
		return this._fc.getLength();
	}
    d.getServer = function()
    {
        return this._server;
    }
	d.leave = function(df)
	{
        console.log('leave('+df+')');
		var o = this._$canvas.parent().offset();
		var self = this;
		
		var $jf = $('<div/>').jFragment({
			'fragment':df.f(),
			'containment':'parent',
			'scroll':false,
			'color':df._fs.fill,
            'distance':0,
			
		}).on('dragstop', function(event, ui) {
                console.log('dragstop called');
				$jf.remove();
                //Tell the server to remove the fragment
                if(df._cf!=undefined)
                {
                   self._server.rmFrag(df._cf.id);
                }
		}).appendTo(this._$canvas.parent());

        
        $jf.data('cf', df._cf);
        
		var l = stage.mouseX - 0.5 * $jf.outerWidth();
		var t = stage.mouseY - 0.5 * $jf.height();

		$jf.css({'position':'absolute', 'left':l, 'top':t,});
		
		this._fc.rm(df);
		
		var jev = $.Event('mousedown', {
            'which':1,
            'pageX':(o.left+stage.mouseX),
            'pageY':(o.top+stage.mouseY),
        });
		console.log('firing mousedown event');
		$jf.trigger(jev);
    }
	
	d.join = function($jf)
	{	
		var f = $jf.jFragment('getFragment')
        var c = $jf.jFragment('getColor');
        var cf = $jf.data('cf');

        //add the fragment into the construct, with dragging
		var df = new DisplayFragment(f,cf);

		df._fs.fill = c;		
		this._fc.add(df);
	}
	
	d._initInfo = function()
	{
		var self = this;
		this._$info = $('.fragment-info');
		
		//setup the buttons
		this._$info.find('#fragment_remove')
			.button({label: 'Remove', icons: {primary:'ui-icon-trash',},})
			.click(function() {self._hideInfo()});

		this._$info.find('#fragment_clip')
			.button({label: 'Clipping', icons: {primary:'ui-icon-scissors'}, disabled:true,})
			.click(function() {});
	}

	d.showInfo = function(df)
	{
		//set the title and description
		this._$info.find('#fragment_name').text(df.f().getName());
		this._$info.find('#fragment_desc').text(df.f().getDesc());
		
		var loc = ra2xy(df.getRadius(), df.getMid());
		loc = this._fc.localToGlobal(loc.x, loc.y);
		
		//set position
		this._$info.css({position: 'absolute', zindex:100, left:loc.x - 33, top:loc.y - (this._$info.outerHeight() + 14),});
		this._$info.css({display:'block',});
		
		//set remove callback
		var self = this;
		this._$info.find('#fragment_remove')
			.unbind('click')
			.click(function() {self.hideInfo();self._fc.rm(df);});
		
		//binding to click means it gets triggered immediately
		this._$canvas.mousedown(function() {self.hideInfo();});
	}

	d.hideInfo = function()
	{
		this._$info.css({display:'none',});
		this._$canvas.unbind('mousedown', this.hideInfo);
	}
	
	
	d._calcSize = function()
	{
		var $c = this._$canvas;
		this._width = $c.width();
		this._height = (10.0 / 16.0) * this._width;
		$c.height(this._height);
		$c.prop('width', this._width);
		$c.prop('height', this._height);
		
		this._child.x = this._width / 2.0;
		this._child.y = this._height / 2.0;
		
		this._tlen.y = 20;
		
		this._server.x = 10;
		this._server.y = this._height - 30;
		
		var radius = Math.min(this._width,this._height) * 0.30; //the base radius of the plasmid
		
		F.radii[Area.CW] = radius + 0.5*F.width;
		F.radii[Area.CCW]= radius - 0.5*F.width;
		
		F.dradii[Area.CW] = F.radii[Area.CW] + 1.6*F.width;
		F.dradii[Area.CCW]= F.radii[Area.CCW] - 1.6*F.width;
			
		F.ldelta[Area.CW] = + 15;
		F.ldelta[Area.CCW]= - 25;
		
		F.maxRadii[Area.CCW] = radius;
		F.maxRadii[Area.CW] = radius + 4.5 * F.width;
		
		F.joinRadius = radius + 4.0 * F.width;
	}
	
	d._gotInfo = function(con)
	{
        //replace with a test construct
        //con = libDesigner.getTestConstruct();

		this.setName(con.name);
	
		this.setLength(con.length);

        var dfs = new Array();
        for(i in con.cfs)
            {
                dfs.push(new DisplayFragment(con.cfs[i].f, con.cfs[i]));
            }
		
		this._fc.addMulti(dfs);
		
		stage.update();
	}
	
	d.toString = function() {return '[Designer '+this._width+'x'+this._height+' (cid='+this._cid+') ]';};
