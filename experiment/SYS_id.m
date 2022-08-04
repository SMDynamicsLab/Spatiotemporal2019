function [ os ] = SYS_id(  )

arq=computer('arch');

if isempty( strfind(arq,'win') )==0
    os=1;
elseif isempty( strfind(arq,'gln') )==0
    os=2;
elseif isempty( strfind(arq,'mac') )==0
    os=3;
end

end
    