/*
 * libDesigner api
 *
 *   depends on ajax.js & libfrag.js
 *
 */

var libDesigner = new function()
{
    this.getConstructByID = function(cid, _suc)
    {
        AJAX.post({
            url: '/gibthon/api/' + cid + '/getInfo/', 
            success: function(data)
            {
                _suc(new Construct(data));
            }
        });
    }
    this.getTestConstruct = function()
    {
        var test_data = {
            id: 1,
            name: 'Test Construct',
            desc: 'This construct is a test',
            length: 10079,
            modified: 'blah',
            fs: [
                {
                "origin": "BioBrick", 
                "name": "pSB1C3", 
                "refs": [], 
                "annots": {}, 
                "length": 2070, 
                "id": 1, 
                "desc": "High copy BioBrick assembly plasmid"
            },
            {
                "origin": "BioBrick", 
                "name": "BBa_K325219", 
                "refs": [], 
                "annots": {}, 
                "length": 3841, 
                "id": 2, 
                "desc": "Red Firefly Luciferase and LRE"+
                    "(under pBAD)L. Cruciata(E. coli optimised)"
            },
            {
                "origin": "Nucleotide Database", 
                "name": "HD065425", 
                "refs": [], 
                "annots": {}, 
                "length": 4168, 
                "id": 5, 
                "desc": "Sequence 52 from Patent WO2010070295."
            },
            ],
            cfs: [
                {
                id: 1,
                direction: 1,
                s_offset: 0,
                s_feat: -1,
                e_offset: 0,
                e_feat: -1,
                order: 0,
            },
            {
                id: 2,
                direction: -1,
                s_offset: 0,
                s_feat: -1,
                e_offset: 0,
                e_feat: -1,
                order: 1,
            },
            {
                id: 3,
                direction: 1,
                s_offset: 0,
                s_feat: -1,
                e_offset: 0,
                e_feat: -1,
                order: 2,
            },
            ],
        };

        return new Construct(test_data);
    };

};

function Construct(data)
{
    this.id = data.id;
    this.name = data.name;
    this.desc = data.desc;
    this.length = data.length;
    this.cfs = new Array();
    this.fs = new Array();
    this.modified = data.created;

    for(var i = 0; i < data.cfs.length; i=i+1)
    {
        var f = new Fragment(data.fs[i]);
        this.fs.push(f);
        this.cfs.push(new ConstructFragment(data.cfs[i], f));
    }

    this.addFragment = function(f, position, direction, _suc)
    {
        console.log('AddingFragment at ' + position);
        AJAX.post({
            url: '/gibthon/api/' + this.id + '/addFragment/', 
            data: {'fid': f.getID(), 'pos': position, 'dir':direction,}, 
            success: function(cf) {
                if($.isFunction(_suc)) _suc(new ConstructFragment(cf, f));
            },
            error: function(jqXHR, textStatus, errorThrown)
            {
                console.log('Could not addFragment: ' + textStatus);
            },
        });
    }

    this.rmFragment = function(cfid, _suc)
    {
        AJAX.post({
            url: '/gibthon/api/' + this.id + '/rmFragment/', 
            data: {'cfid': cfid,}, 
            success: function() {if(_suc!=undefined) _suc();},
        });
    }

    this.reorder = function(cfids, dirs, _suc)
    {
        AJAX.post({
            url: '/gibthon/api/' + this.id + '/saveOrder/', 
            data: {'cfid[]':cfids, 'direction[]':dirs,},
            success: function() {if(_suc!=undefined) _suc();},
        });
    }

    this.saveMeta = function(name, desc)
    {
        if(name)//can't have an empty string for name
            this.name = name;
        if(desc!=undefined)
            this.desc = desc;
        AJAX.post({
            url: '/gibthon/api/' + this.id + '/saveMeta/', 
            data: {'name': this.name, 'desc': this.desc,}, 
        });
    }



}

function ConstructFragment(d, f)
{		
    this.id = d.id;
    this.strand = d.direction;
    this.s_offset = d.s_offset;
    this.s_feat = d.s_feat;
    this.e_offset = d.e_offset;
    this.e_feat = d.e_feat;
    this.order = d.order;
    this.f = f;

    this.startPos = function()
    {
        /*if(this.id == undefined) return 0;
        //console.log('startPos: s_feat: '+this.s_feat+' s_offset: '+this.s_offset);
        if(this.s_feat > 0)
        {
        //var sf = f.getFeatById(this.s_feat);
        if(sf != null)
        return this.s_offset + sf.start;
        }*/
        return this.s_offset;
    };
    this.endPos = function()
    {
        /*if(this.id == undefined) return this.f.getLength();
          if(this.e_feat > 0)
          {
          var ef = this.f.getFeatById(this.e_feat);
          if(ef != null)
          return this.e_offset + ef.end;
          }*/
        return this.f.getLength() - this.s_offset;
    };
    this.getLength = function()
    {
        return Math.abs(this.endPos() - this.startPos());
    };
    this.length = function()
    {
        console.warn('ConstructFragment.length depreciated');
        return this.getLength();
    }
    this.toString = function() {return '[ConstructFragment (id='+this.id+') ]';};
};

/*
 * * jQuery widget for previewing the designed fragment
 * */
$.widget('ui.constructPreview', {
    options: {
        primers: {},
    },
    _create: function() {
        var self = this;
        var el = $(this.element[0]);
        this.el = el;
        this.pview = el.find('.primer-view');
        this.arrow = this.pview.find('.arrow');
        this.pedit = el.find('.primer-edit');
        this.boxplot = el.find('.boxplot');
        this.type = 'fwd';

        var css = $("<link>");
        css.attr({
            rel: "stylesheet",
            type: "text/css",
            href: "/static/css/construct_preview.css"
        });
        $("head").append(css);

        el.find('#overview .primer').click( function(){
            if($(this).hasClass('sl'))
            {
                el.find('.primer.sl').removeClass('sl');
                self.pview.slideUp();
                return;
            }
            el.find('.primer.sl').removeClass('sl');
            $(this).addClass('sl');
            self.show($(this));
        });

        this.boxplot.dialog({
            autoOpen: false,
            title: 'Boxplot',
            resizable: false,
            modal: true,
            width: 700,
            height: 850,
        });

        el.find('button.boxplot-btn').button({
            icons: {primary: 'ui-icon-image',},
        }).click( function() {
            self.boxplot.html('')
                .dialog('open')
                .load('primers/' + $(this).attr('pid') + '/boxplot/');
        });

        $(document).on('click', '.join .seq > .base', function() {
            $(this)
            .addClass('sl')
            .siblings()
            .removeClass('sl');
            self._set_primer_pos();
        });
    },
    _init: function() {
        this.el.find('#overview .fragment:not(.bk-fragment)').each( function(i) {
            $(this).find('.seq').css({
                'background-color': libFrag.getNextColor(),
            });
        });

    },
    show: function($p)
    {
        this.arrow.css({'left': this._get_center($p), });
        $p.hasClass('fwd-primer') ? this.type='fwd' : this.type='rev';
        var p = this.options.primers[$p.attr('id')];
        this.pedit.find('#pname').text(p.name);
        this.pedit.find('#pseq').text(p.seq);
        
        if(this.type == 'fwd')
        {
            this._set_fwd();
            this.pedit.find('#left_fragment > .fname').text(
                this.el.find('#fragment-' + p.flap.cfid + ' > .fname').text());
            this.pedit.find('#right_fragment > .fname').text(
                this.el.find('#fragment-' + p.stick.cfid + ' > .fname').text());
            this._fwd_seq(p.flap.context, p.stick.context);
            this.pedit
                .find('#left_fragment .base[offset="' + (p.flap.length-1) + '"]')
                .addClass('sl');
            this.pedit
                .find('#right_fragment .base[offset="' + (p.stick.length-1) + '"]')
                .addClass('sl');
        }
        else
        {
            this._set_rev();
            this.pedit.find('#left_fragment > .fname').text(
                this.el.find('#fragment-' + p.stick.cfid + ' > .fname').text());
            this.pedit.find('#right_fragment > .fname').text(
                this.el.find('#fragment-' + p.flap.cfid + ' > .fname').text());
            this._rev_seq(p.flap.context, p.stick.context);
            this.pedit
                .find('#left_fragment .base[offset=' + (p.stick.length-1) + ']')
                .addClass('sl');
            this.pedit
                .find('#right_fragment .base[offset=' + (p.flap.length-1) + ']')
                .addClass('sl');
        }

        this.pedit.find('.prule-full > .text').html(this._ruler_full(p));
        this.pedit.find('.prule-stick > .text').html(this._ruler_stick(p));

        var w = this.pedit.find('.primer-warnings > ul').html('');
        if(p.warnings.length == 0)
            w.parent().hide();
        else
        {
            w.parent().show();
            for(var i = 0; i < p.warnings.length; i=i+1)
            {
                $('<li>' + p.warnings[i] + '</li>').appendTo(w);
            }
        }

        this.pedit.find('.boxplot-btn').attr('pid', p.id);
        this._set_colors($p);
        this.pview.slideDown();
        this._set_primer_pos();
    },
    _set_primer_pos: function()
    {
        var l = $('#left_fragment .base.sl');
        var r = $('#right_fragment .base.sl');
        p = this.options.primers[this.pedit.find('#pname').text()];
        if(this.type == 'fwd')
        {
            this.pview.find('.fwd-primer').css({
                'left': l.position().left - l.parent().width(),
                'width': r.offset().left - l.offset().left + r.outerWidth(),
            }).find('.prule-stick').css({
                'left': l.parent().width() - l.position().left,
            });
        }
        else
        {
            this.pview.find('.rev-primer').css({
                'left': l.position().left,
                'width': r.offset().left - l.offset().left + r.outerWidth(),
            }).find('.prule-stick').css({
                'right': r.position().left + r.width(),
            });
        }
    },
    _get_center: function(f)
    {
        return f.offset().left - this.el.offset().left + 15 + 0.5*f.outerWidth();
    },
    _set_colors: function(p)
    {
        var l, r;
        if(p.hasClass('fwd-primer'))
        {
            r = p.siblings('.seq').css('background-color');
            var prev = p.parent().prev('.fragment:not(.bk-fragment)');
            if(prev.length == 0)
                prev = p.parent().siblings('.fragment:not(.bk-fragment)')
                    .last();
            l = prev.find('.seq').css('background-color');
        }
        else
        {
            l = p.siblings('.seq').css('background-color');
            var next = p.parent().next('.fragment:not(.bk-fragment)');
            if(next.length == 0)
                next = p.parent().siblings('.fragment:not(.bk-fragment)')
                    .first();
            r = next.find('.seq').css('background-color');
        }

        this.pview.find('#left_fragment .seq').css('background-color', l);
        this.pview.find('#right_fragment .seq').css('background-color', r);
    },
    _set_fwd: function()
    {
        this.pedit.find('.rev-primer').hide();
        this.pedit.find('.fwd-primer').show();
    },
    _set_rev: function()
    {
        this.pedit.find('.fwd-primer').hide();
        this.pedit.find('.rev-primer').show();
    },
    _fwd_seq: function(flap, stick)
    {
        this.pedit.find('.seq').removeClass('has_seq').html('');
        this.pedit.find('.fwd_seq').addClass('has_seq');

        var s = this.pedit.find('#right_fragment .fwd_seq');
        var f = this.pedit.find('#left_fragment .fwd_seq');

        for(var i = 0; i < stick.length; i = i+1)
            $("<div>").addClass('base').attr('offset', i)
                .text(stick[i]).appendTo(s);
        for(var i = 0; i < flap.length; i=i+1)
            $("<div>").addClass('base').attr('offset', i)
                .text(flap[flap.length -1 - i]).appendTo(f);

    },
    _rev_seq: function(flap, stick)
    {
        this.pedit.find('.seq').removeClass('has_seq').html('');
        this.pedit.find('.rev_seq').addClass('has_seq');

        var s = this.pedit.find('#left_fragment .rev_seq');
        var f = this.pedit.find('#right_fragment .rev_seq');

        for(var i = 0; i < stick.length; i = i+1)
            $("<div>").addClass('base').attr('offset', i)
                .text(stick[i]).appendTo(s);
        for(var i = 0; i < flap.length; i=i+1)
            $("<div>").addClass('base').attr('offset', i)
                .text(flap[flap.length - 1 -i]).appendTo(f);
    },
    _ruler_full: function(p)
    {
        return p.length + "bp (" + p.tm + "&deg;C)";
    },
    _ruler_stick: function(p)
    {
        return p.stick.length + "bp (" + p.stick.tm + "&deg;C)";
    },
});

