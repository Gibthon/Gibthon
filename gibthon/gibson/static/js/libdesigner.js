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

    },
    _create: function() {
        var self = this;
        var el = $(this.element[0]);
        this.el = el;
        this.pview = el.find('.primer-view');
        this.arrow = this.pview.find('.arrow');
        this.pedit = el.find('.primer-edit');
        this.boxplot = el.find('.boxplot');

        var css = $("<link>");
        css.attr({
            rel: "stylesheet",
            type: "text/css",
            href: "/static/css/construct_preview.css"
        });
        $("head").append(css);

        el.find('.primer').click( function(){
            if($(this).hasClass('sl'))
            {
                el.find('.primer.sl').removeClass('sl');
                self.pview.slideUp();
                return;
            }
            el.find('.primer.sl').removeClass('sl');
            $(this).addClass('sl');
            self.arrow.css({'left': self._get_center($(this)), });

            var type;
            $(this).hasClass('fwd-primer') ? type = 'fwd' : type = 'rev';
            
            self.pedit.html('');
            console.log('$(this).find(".primer-detail").length = ' + 
                       $(this).find('.primer-detail').length);
            $(this).find('.primer-detail')
                .children().clone(true).appendTo(self.pedit);

            self.pview.slideDown();
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

        $(document).on('click', '.join > .seq > .base', function() {
            $(this)
            .addClass('sl')
            .siblings()
            .removeClass('sl');
        });

        $(document).on('click', '#right_seq > .base', function() {
            var t = $(this);
            self.pview.find('.primer').css({
                'right': t.offsetParent().width() - 
                    t.position().left - t.outerWidth(),
            });
        });
        $(document).on('click', '#left_seq > .base', function() {
            self.pview.find('.primer').css({
                'left': $(this).position().left,
            });
        });



    },
    _init: function() {
        this.el.find('.fragment:not(.bk-fragment)').each( function(i) {
            $(this).find('.seq').css({
                'background-color': libFrag.getNextColor(),
            });
        });

    },
    loadJson: function(json)
    {
        var p = this.pview.find('.join > .primer').prop('id', 'primer-'+json.id);
        var l, r, l_seq, r_seq;
        if(json.type == 'fwd')
        {
            l = json.flap_len;
            l_seq = json.flap_seq;
            r = json.stick_len;
            r_seq = json.stick_seq;            
            p.removeClass('rev-primer').addClass('fwd-primer');
        }
        else
        {
            l = json.stick_len;
            l_seq = json.stick_seq;
            r = json.flap_len;
            r_seq = json.flap_seq;
            p.removeClass('fwd-primer').addClass('rev-primer');
        }

        this._fill(this.pview.find('#left_seq'), l_seq);
        this._fill(this.pview.find('#right_seq'), r_seq);

        var lb = this.pview.find('.join > #left_seq > [offset='+l+']').addClass('sl');
        var rb = this.pview.find('.join > #right_seq > [offset='+r+']').addClass('sl');

        console.log('lb.position().left = ' + lb.position().left);
        console.log('rb.position().left = ' + rb.position().left);

        this.pview.find('.primer').css({'left': lb.position().left,});
        this.pview.find('.primer').css({
            'right': rb.offsetParent().width() - rb.position().left - rb.outerWidth(),
        });    

    },
    toJson: function()
    {
        var l = this.pview.find('.join > #left_seq  > .sl').attr('offset');
        var r = this.pview.find('.join > #right_seq > .sl').attr('offset');

        var s;
        var f;

        var type;

        if($('.join > .primer').hasClass('fwd-primer'))
        {
            f = l;
            s = r;
            type = 'fwd';
        }
        else
        {
            f = r;
            s = l;
            type = 'rev';
        }
        return {'id': $('.join > .primer').prop('id').split('-')[1],
            'stick_len': s,
            'flap_len': f,
            'type': type,};
    },
    _get_center: function(f)
    {
        return f.offset().left - this.el.offset().left + f.outerWidth();
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

        this.pview.find('#left_seq').css('background-color', l);
        this.pview.find('#right_seq').css('background-color', r);
    },
    _fill: function($d, seq)
    {
        $d.html('');
        if($d.prop('id') == 'left_seq')
        {
            seq = seq.split('').reverse().join('');
        }
        for(var b = 0; b < seq.length; b=b+1)
        {
            $('<div>')
            .addClass('base')
            .text(seq[b])
            .attr('offset', b+1)
            .appendTo($d);
        }

        $('<div>')
            .css({'clear':'both',})
            .appendTo($d);
    },
});

