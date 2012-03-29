/*
 * designer.js, Haydn King - hjking734@gmail.com - 2012
 * 
 * */

//Global things

	//the stage
	var stage;


	/**
	 * Areas - represent the possible areas that a fragment could be in
	 * 
	**/
	var Area = {
		CCW: 0,
		CW: 1,
		RM: 2,
	};

	var AreaString = ['ccw', 'cw', 'rm',];

	var a2s = function(a) {return AreaString[a];};

//angle Functions
	var _D_PER_R = 360.0 / (2 * Math.PI);
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
		F.radii = [0,0,0,];
		/**
		 * Stores the radii for each area when a fragment is being dragged
		 * @property dradii
		 * @public
		 **/
		F.dradii = [0,0,0,];
		/**
		 * Stores the label radii for each area
		 * @property lradii
		 * @public
		 **/
		F.lradii = [0,0,0,];
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
		 * The maximum radius in the construct
		 * @property maxRadius
		 * @public
		 **/
		F.maxRadius = 0;
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
		console.log('fs.draw(...)');
		p = ['width', 'angle', 'rotation', 'radius'];
		for(i in p)
			console.log('	fs.'+p[i]+' = ' + this[p[i]]);
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
		if(_update == undefined) _update = true;
		
		this.text = _text;
		this._redraw();
		//Cache
		
		if(_update)
			stage.update();
		return this;
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
	fl.setRadius = function(_radius, _update)
	{	
		if(_update == undefined) _update = true;
		
		this._radius = _radius;
		this._realign();
		
		if(_update)
			stage.update();
		return this;
	}
	
	/**
	 * Get the rotation
	 * @method getRadius
	 **/
	fl.getRotation = function()
	{
		return deg2rad(this.rotation);
	};
	
	/**
	 * Set the rotation
	 * @method setRotation
	 * @param {float} _rotation The new rotation (rads)
	 * @param {boolean} _update Whether to update the stage
	 * @default true
	 **/
	fl.setRotation = function(_rotation, _update)
	{	
		if(_update == undefined) _update = true;
		
		this.rotation = r2d(_rotation);
		
		if(_update)
			stage.update();
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
		//for each letter
		for(var i = 0; i < this.getNumChildren(); i = i + 1)
		{
			var l = this.getChildAt(i);
			//reset the center
			l.x = 0;
			//set the y-offset accordingly
			if(this._outward)
			{
				l.y = this.radius;
				l.regY = -this._radius;
			}
			else
			{
				l.y = -this.radius;
				l.regY = this._radius;
			}
				
			//set the rotation
			l.rotation = r2d(r);
			if(!this.outward)
				l.rotation = -r2d(r);
			
			r = r + l.getMeasuredWidth() / this.radius;
		}
		return this;
	};
	
/**
* 
* @class DisplayFragment
* @extends Container
* @constructor
* @param {Fragment} _f The Fragment which this represents
* @param {ConstructFragment} _df The ConstructFragment - represents it's position in this construct
* @default null
* @param {Area} _area The initial area
* @default Area.RM
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
		this._f = _f;
		this._cf = _cf;
		this._area = Area.RM;
		if(_cf != undefined)
		{
			if(_cf.strand = 1)
				this._area = Area.CW;
			if(_cf.strand = -1)
				this._area = Area.CCW;
		}

		this._fs = new FragmentShape(0,0,F.radii[this._area]);
		this._fl = new FragmentLabel(this._f.name, F.lradii[this._area], this._fs.angle / 2.0);
		
		this.addChild(this._fs);
		this.addChild(this._fl);
		
		if(this._area == Area.RM)
			this._fl.hide();
		
		/* For now, don't support immediate dragging
		if(_df == undefined)
		{
			this._df = null;
			//we are dragging
			this._drag = true;
		}
		else
		{
			this.drag = false;
			//not dragging
		}
		* */
		this._drag = false;		
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
		if(_a == RM)
		{
			this._fs.linear = true;
			this._fl.hide();
		}
		else
		{
			this.x = 0; this.y = 0;
			
			this._fs.linear = false;
			if(this._drag)
			{
				this._fs.radius = F.radii[_a];
				this._fl.hide();
			}
			else
			{
				this._fs.radius = f.dradii[_a];
				this._fl.show();
			}
			if(_a == CW)
				this._fl.clockwise = 1.0;
			else
				this._fl.clockwise = -1.0;
		}
	}
	 
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
		if(this._cf)
			return this._cf.length();
		return this._f.length;
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
	
	/**
	 * Get the angle (rads) at which the fragment starts
	 * @method getStart
	 **/
	df.getStart = function()
	{
		return d2r(this._df.rotation);
	}
	
	/**
	 * Get the angle (rads) at the middle
	 * @method getMid
	 **/
	df.getMid = function()
	{
		return d2r(this._cf.rotation) + this._fs.angle / 2.0;
	}
	
	/**
	 * Get the angle (rads) at which the fragment ends
	 * @method getEnd
	 **/
	df.getEnd = function()
	{
		return d2r(this._df.rotation) + this._fs.angle;
	}
	
	/**
	 * Get the angle (rads) of the fragment (i.e. length)
	 * @method getAngle
	 **/
	df.getAngle = function()
	{
		return this._fs.angle;
	}
	
	/**
	 * Animate the properties of the framgment
	 * @method animate
	 * @param {Object} t The target properties
	 * Supports start (rads), angle(rads), radius, alpha, clockwise
	 * @param {function} cb A callback for when the animation completes 
	 **/
	df.animate = function(t, cb)
	{
		var self = this;
		//make sure any rotation goes the fastest route
		if(t.rotation != undefined)
		{
			var r = bound_degs(r2d(t.rotation) - this._fs.rotation); //r is the distance and direction of shortest rotation
			t.rotation = this._fs.rotation + r;
		}
		
		var done = function() {
			this._fs.rotation = bound_angle(self.rotation);
			if((this._area != RM) && (!this._drag))
				this._fl.show();
		};
		
		this._fl.hide();
		var tween = Tween.get(this._fs)
		 .call(setAnim) //enable global animation
		 .to(t, 250, Ease.quartOut)
		 .call(done)
		 .call(clearAnim); //disable global animation
		if(cb != undefined) tween.call(cb);
	}
	
	/**
	 * Set the properties of the framgment
	 * @method set
	 * @param {Object} t The target properties
	 * Supports start (rads), angle(rads), radius, alpha, clockwise
	 **/
	df.set = function(t)
	{
		console.log('  df.set({angle: '+t.angle+', rotation: '+t.rotation+'})');
		//set rotation
		if(t.rotation != undefined)
		{
			this._fs.rotation = r2d(bound_rads(t.rotation));
		}
		//set angle
		if(t.angle != undefined)
		{
			this._fs.angle = bound_rads(t.angle);
		}
		//set radius
		if(t.radius != undefined)
		{
			this._fs.radius = t.radius;
		}
		//set alpha
		if(t.alpha != undefined)
		{
			this._fs.alpha = t.alpha;
		}
		//set clockwise
		if(t.clockwise != undefined)
		{
			this._fs.clockwise = t.clockwise;
		}		
	}
	
	/**
	 * Called just before draw() - update label position
	 * @method tick
	 **/
	df.tick = function()
	{
		this._fl.setRotation(this.getMid(), false);
	}
	
	df.toString = function()
	{
		return '[DisplayFragment (f.name = "'+this._f.name+'")]'; 
	}
	
//Mouse Events
	/**
	 * Handle a mouse Over
	 * @method onMouseOver
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onMouseOver = function(ev)
	{
		console.log('"'+this.f.name+'".onMouseOver('+ev.stageX+','+ev.stageY+')');
		this.alpha = 0.8;
		stage.update();
	}
	
	/**
	 * Handle a mouse Out
	 * @method onMouseOut
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onMouseOut = function(ev)
	{
		console.log('"'+this.f.name+'".onMouseOut('+ev.stageX+','+ev.stageY+')');
		this.alpha = 1.0;
		stage.update();
	}
	
	//timeout for a drag
	var dt = null;
	
	/**
	 * Handle a mouse click
	 * @method onClick
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onClick = function(ev)
	{
		console.log('"'+this.f.name+'".onMouseClick('+ev.stageX+','+ev.stageY+')');
		//if we were about to start dragging
		if(dt != null)
		{
			clearTimeout(dt);
			dt = null;
			//do on-click stuff
			
			return;
		}
		//onDrop should have already been called
	}
	
	/**
	 * Handle a mouse press
	 * @method onPress
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onPress = function(ev)
	{
		console.log('"'+this.f.name+'".onPress('+ev.stageX+','+ev.stageY+')');
		
		//set a timeout to handle dragging if the mouse isn't released
		var self = this;
		dt = setTimeout( function() {
			dt = null;
			self._drag = true;
			self.animate({radius: F.dradii[this.area],});
			
			//hookup the events
			ev.onMouseMove = self.onDrag;
			ev.onMouseUp = self.onDrop;
		}, 250);
	}
	
	/**
	 * Handle a mouse Drag
	 * @method onDrag
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onDrag = function(ev)
	{
		console.log('"'+this.f.name+'".onDrag('+ev.stageX+','+ev.stageY+')');
		var p = this._get_mev(ev);
		this._rotate(p.a);
		this.parent.SortOne(this);
	};
	
	/**
	 * Handle a mouse drop
	 * @method onDrop
	 * @param {MouseEvent} ev The mouse event
	 **/
	df.onDrop = function(ev)
	{
		console.log('"'+this.f.name+'".onDrop('+ev.stageX+','+ev.stageY+')');
		this._drag = false;
		this.parent.onDrop(this);
	}
	

//Private Methods
	/**
	 * Set the rotation, given the mouse is at angle a (rads)
	 * @method _rotate
	 * @private
	 * @param {angle} a The angle of the mouse
	 **/
	df._rotate = function(a)
	{
		self._fs.rotation = r2d(a - this._mouse_offset);
	}
	
	/**
	 * Get the local radial and euclidian mouse position
	 * @method _get_mev
	 * @private
	 * @param {MouseEvent} ev The MouseEvent
	 **/
	df._get_mev = function(ev)
	{
		var p = this.globalToLocal(ex.stageX, ev.stageY);
		var rp = xy2ra(p.x,p.y);
		return { x: p.x, y:p.y, r: rp.r, a: rp.a,};
	}


/**
* Manitains the position of the fragments while they're in the construct
* @class FragmentContainer
* @extends Container
* @constructor
* @param {Server} server A server object for saving changes
**/
var FragmentContainer = function(server)
{
	this.initialize(server);
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
	fc.initialize = function(server) 
	{
		this.Container_initialize();
		this._server = server;
	};
	
//public methods
	/**
	 * Add a new DisplayFragment
	 * @method addFragAt
	 * @param {DisplayFragment} df The DisplayFragment to add
	 * @param {int} pos The position to add the DisplayFragment
	 * @default df.cf.order
	 **/
	fc.addFragAt = function(df, pos)
	{
		if(pos == undefined)
			pos = df.cf.order;
		else
			pos = this._bound(pos);
				
		this.addChildAt(df, pos);
		for(var i = pos; i < this.getNumChildren(); i = i+1)
		{
			var c = this.getChildAt(i);
			if(c.df != undefined)
				c.df.order = i;
		}
		
		this._update_length();
		
		this._update_layout( 
			pos+1,
			this.getFragAt(pos+1).getStart() + 2*Math.Pi*cf.getLength()/this._eff_length(),
			true
		);
		
		return this;
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
		
		for(i in this.children)
		{
			this.children[i].getCf().order = i;
		}
		
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
		return this;
	}
	
	/**
	 * Notify the container that a drop occured
	 * @method drop
	 * @param {DisplayFragment} df
	 **/
	fc.onDrop = function(df)
	{
		this._update_layout();
		return this;
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
			startFrag = this._bound(startFrag);
			
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
			
			if(!df.drag)
				t.rotation = r;
			
			//calculate the angle
			var l = df.getLength();
			if(l < this._lmin) l = this._lmin;
			
			var a = 2 * Math.PI *  l / this._eff_length;
			
			console.log('calculate angle: var a = '+a+' = 2 * Math.PI * '+df.getLength()+' / '+this._eff_length+';');
			
			t.angle = a;
			console.log('t.angle = '+t.angle+' = a = '+a+';');
			r = r + a;
			console.log('push target = {angle: '+t.angle+', rotation: '+t.rotation+'}');
			targets.push(t);
		}
		
		console.log('fc._updateLayout(startFrag = '+startFrag+',startAngle = '+startAngle+',animate = '+animate+');');
		console.log('  result:');
		for(t in targets)
		{
			var f = this.getFragAt(t + startFrag);
			console.log(this._bound(t+startFrag) + ' ' + targets[t].rotation + ' -> ' + (targets[t].rotation + targets[t].angle));
		}
			
		//apply the targets all at once
		if(animate)
		{
			for(t in targets)
				this.getFragAt(t + startFrag).animate(targets[t]);
		}
		else
		{
			for(t in targets)
			{
				console.log('this.getFragAt('+ t +').set(targets['+t+'] = {rotation: '+targets[t].rotation+', angle: '+targets[t].angle + '});');
				this.getFragAt(t + startFrag).set(targets[t]);
			}
		}
		
		return this;
	}
	
	/**
	 * sort the fragment s, assuming all others are in order
	 * @method _sortOne
	 * @param {int} s
	 * @public
	 **/
	fc.sortOne = function(s)
	{
//TODO: implement bi-directional bubble...
		this._sortAll();
		return this;
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
		
		//find the minimum length
		this._lmin = Math.ceil(this._length * F.minAngle / (2 * Math.PI));
		
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
		while(i < 0) i = i + this.getNumChildren;
		return i % this.getNumChildren();
	}

	
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


	s._container_init = s.initialize;

	s.initialize = function()
	{
		this._text.x = 0;
		this._sm_text.x = 10;
		this._sm_text.y = -this._text.getMeasuredLineHeight() - 2;
		this._sm_text.visible = false;
		this._text.visible = false;
		
		this.addChild(this._text, this._sm_text);
	}

	//Public Functions

	s.getInfo = function(cb, cid)
	{
		var self = this;
		cd_get_info(function(data)
			{
				//data = {name, desc, Array[fragment] fs, Array[construct fragment] cfs}
				if(data.fs.length != data.cfs.length)
				{
					self._handle_error('getting info', "Couldn't parse response", "number of fragment != number of constructFragments");
					return;
				}
				var dfs = new Array();
				for(var i=0; i < data.fs.length; i=i+1)
				{
					dfs.push( new DisplayFragment(data.fs[i], new ConstructFragment(data.cfs[i], data.fs[i])) );
				}
				
				cb(data.name, data.desc, data.length, dfs);
				
			}, cid, this._handle_error);
		return this;
	}

	//Private functions
	
	s._handle_error = function(desc, textStatus, errorThrown)
	{
		this._text.text = "Sorry, there was an error while " + desc + ", please try reloading.";
		this._sm_text.text = "status: '" + textStatus + "' error: '" + errorThrown + '"';
		this._text.color = UI.RED;
		this._sm_text.color = UI.BLACK;
		this._sm_text.visible = true;
		this.visible = true;
		stage.update();
	}
	
	s._handle_success = function(cb, data)
	{
		this.visible = false;
		this._sm_text.visible = false;
		cb(data);
	}

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
		this.$canvas = $canvas;
		this._cid = cid;
		
		this._child = new Container();
		this._tname = new Text('', FONT.LG, COL.BLACK);
		this._tlen = new Text('', FONT.M, COL.BLACK);
		this._tname.textAlign = 'center';
		this._tname.maxWidth = 1000;
		this._tlen.textAlign = 'center';
		this._tlen.maxWidth = 1000;
		
		this._fc = new FragmentContainer();
		
		this._server = new Server();
		
		this._child.addChild(this._tname, this._tlen, this._fc);
		this.addChild(this._child, this._server);
		
		this._calc_size();
		
		this.setName('My Name (' + this._cid + ')');
		this.setLength(150);
		
		//setup the canvas
		stage = new Stage(document.getElementById('cdesigner')); //$canvas.get());
		stage.addChild(this);
		
		Ticker.setFPS(25);
		Ticker.addListener(this);
		
		var self = this;
		this._server.getInfo(function(a,b,c,d) {self._gotInfo(a,b,c,d)}, this._cid);
		
		stage.update();
	}
	
	//Public methods
	d.tick = function()
	{
		// update.
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
	
	//private methods
	d._calc_size = function()
	{
		var $c = this.$canvas;
		this._width = $c.width();
		this._height = (7.0 / 16.0) * this._width;
		$c.height(this._height);
		$c.prop('width', this._width);
		$c.prop('height', this._height);
		
		this._child.x = this._width / 2.0;
		this._child.y = this._height / 2.0;
		
		this._tlen.y = 20;
		
		this._server.x = 10;
		this._server.y = this._height - 30;
		
		var radius = Math.min(this._width,this._height) * 0.35; //the base radius of the plasmid
		
		F.radii[Area.CW] = radius + F.width;
		F.radii[Area.CCW]= radius - F.width;
		F.radii[Area.RM] = 1000;
		
		F.dradii[Area.CW] = F.radii[Area.CW] + F.width;
		F.dradii[Area.CCW]= F.radii[Area.CCW] - F.width;
		F.dradii[Area.RM] = F.radii[Area.RM];
			
		F.lradii[Area.CW] = radius + 1.5 * F.width;
		F.lradii[Area.CCW]= radius - 1.5 * F.width;
		F.lradii[Area.RM] = 1000;
		
		F.maxRadius = radius + 2.5 * F.width;
	}
	
	d._gotInfo = function(name, desc, length, dfs)
	{
		this.setName(name);
	
		this.setLength(length);
		
		this._fc.addMulti(dfs);
		
		stage.update();
	}
