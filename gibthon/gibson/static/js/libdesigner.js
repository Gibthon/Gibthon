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
        this.cfs.push(new ConstructFragment(cfs[i], f));
    }

   this.addFragment = function(fid, position, direction, _suc)
    {
        AJAX.post({
            url: '/gibthon/api/' + this.id + '/addFragment/', 
            data: {'fid': fid, 'pos': position, 'dir':direction,}, 
            success: function() {if(_suc!=undefined) _suc();},
            error: function(jqXHR, textStatus, errorThrown)
            {
                console.log('Could not addFragment: ' + textStatus);
            },
        });
    }

    this.rmFragment = function(fid, _suc)
    {
        AJAX.post({
            url: '/gibthon/api/' + cid + '/rmFragment/', 
            data: {'cfid': cfid,}, 
            success: function() {if(_suc!=undefined) _suc();},
        });
    }

    this.reorder = function(fids, _suc)
    {
        ajax.post({
            url: '/gibthon/api/' + cid + '/saveOrder/', 
            data: {'d[]': d,},
            success: function() {if(_suc!=undefined) _suc();},
        });
    }

    this.saveInfo = function()
    {
        ajax.post({
            url: '/gibthon/api/' + cid + '/saveMeta/', 
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
		if(this.id == undefined) return 0;
		//console.log('startPos: s_feat: '+this.s_feat+' s_offset: '+this.s_offset);
		if(this.s_feat > 0)
		{
			var sf = f.getFeatById(this.s_feat);
			if(sf != null)
				return this.s_offset + sf.start;
		}
		return this.s_offset;
	};
	this.endPos = function()
	{
		if(this.id == undefined) return f.length;
		if(this.e_feat > 0)
		{
			var ef = f.getFeatById(this.e_feat);
			if(ef != null)
				return this.e_offset + ef.end;
		}
		return f.length - this.s_offset;
	};
	this.length = function()
	{
		return Math.abs(this.endPos() - this.startPos());
	};
	this.toString = function() {return '[ConstructFragment (id='+this.id+') ]';};

 }
