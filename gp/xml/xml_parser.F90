module xml_parser
  use xml_structure
  implicit none

  private
  public tag
  public xml_open,xml_close
  public xml_read,xml_write
  public get_xml_data,add_xml_data
  public tag_open,tag_close,tag_next
  public name,attribute

  
  interface tag_next
     module procedure get_next_,get_next_named
  end interface

  interface add_xml_data
     module procedure add_real_data, add_real_vector_data, add_char_data, add_char_data_vector
  end interface

  interface get_xml_data
     module procedure get_real_data, get_char_data  !,get_real_vector_data
  end interface

  interface attribute
     module procedure get_attribute,set_attribute
  end interface

  interface name
     module procedure get_name,set_name
  end interface
  
! A type for the enduser in which the data structure is hidden.
  
  include "xml_settings.h"

  type(tag),pointer,save     :: xml_documents(:)
  type(tag),pointer,save     :: xml_index
  type(tag_ref),pointer,save :: nodes(:)
  integer,save               :: nodes_index


contains
  
  subroutine xml_init()
    nullify(xml_documents)
  end subroutine


  function tag_open(cur_tag,name,id)
    character(len=*) :: name
    character(len=*),optional :: id
    type(tag),pointer :: cur_tag
    type(tag),pointer :: tag_open

    if (.not.(associated(xml_documents))) cur_tag=xml_open('new')  
    tag_open=>add_child(cur_tag,name)
    if (present(id)) call add_attr(tag_open,'name',id)
  end function

  function tag_close(cur_tag)
    type(tag),pointer :: cur_tag
    type(tag),pointer :: tag_close
    if(associated(cur_tag%parent)) then; tag_close=>cur_tag%parent    
    else; tag_close=>null(); end if
  end function
  
  subroutine set_attribute(cur_tag,name,value)
    type(tag),pointer :: cur_tag
    character(len=*) :: name,value
    call add_attr(cur_tag,name,value)
  end subroutine

  character(len=attr_len) function get_attribute(cur_tag,name)
    type(tag),pointer :: cur_tag
    character(len=*) :: name
    get_attribute=get_attr(cur_tag,name)
  end function

  subroutine set_name(cur_tag,name)
! Changes tagname
    type(tag),pointer :: cur_tag
    character(len=*) :: name
    cur_tag%name=name
  end subroutine

  character(len=name_len) function get_name(cur_tag)
    type(tag),pointer :: cur_tag
    get_name=cur_tag%name
  end function


  subroutine add_real_data(cur_tag,real_data)
    !type(tag),pointer,intent(out) :: cur_tag
    type(tag),pointer :: cur_tag
    real(kind=real_kind),intent(in) :: real_data(:,:)

    allocate(cur_tag%real_data(size(real_data,1),size(real_data,2)))
    cur_tag%real_data=real_data
    call add_attr(cur_tag,'type','real_data')
  end subroutine

  subroutine add_real_vector_data(cur_tag,real_data)
    !type(tag),pointer,intent(out) :: cur_tag
    type(tag),pointer :: cur_tag
    real(kind=real_kind),intent(in) :: real_data(:)
    allocate(cur_tag%real_data(1,size(real_data)))
    cur_tag%real_data(1,:)=real_data
    call add_attr(cur_tag,'type','real_data')
  end subroutine

  subroutine add_char_data(cur_tag,char_data)
    !type(tag),pointer,intent(out) :: cur_tag
    type(tag),pointer :: cur_tag
    character(len=*),intent(in) :: char_data
    call add_data(cur_tag,char_data)
  end subroutine

  subroutine add_char_data_vector(cur_tag,char_data)
    !type(tag),pointer,intent(out) :: cur_tag
    type(tag),pointer :: cur_tag
    character(len=*),intent(in) :: char_data(:)
    integer :: i
    do i=1,size(char_data)
       call add_data(cur_tag,char_data(i))
    end do
  end subroutine


  subroutine get_real_data(tmp_data,cur_tag)
    type(tag),intent(in) :: cur_tag
   ! real(kind=8),pointer,intent(out) :: tmp_data(:,:)
    real(kind=8),pointer :: tmp_data(:,:)

    tmp_data=>null()
    if (associated(cur_tag%real_data)) tmp_data=>cur_tag%real_data       
  end subroutine

!  function get_real_vector_data() result(tmp_data)
!    real(kind=real_kind),pointer :: tmp_data(:)
!
!    tmp_data=>null()
!    if (.not.(associated(nodes))) return
!    if (nodes_index>=size(nodes)) return
!    if (associated(nodes(nodes_index)%node%real_data)) tmp_data=>nodes(nodes_index)%node%real_data(1,:)       
!  end function

  subroutine get_char_data(tmp_data,cur_tag) 
    type(tag),intent(in) :: cur_tag
    !character(len=data_len),pointer,intent(out) :: tmp_data(:)
    character(len=data_len),pointer :: tmp_data(:)
    tmp_data=>null()
    if (associated(cur_tag%data)) tmp_data=>cur_tag%data 
  end subroutine

  integer function get_next_(cur_tag)
    type(tag),pointer :: cur_tag
    integer :: i
    type(tag_ref),pointer :: nodes(:)
    if(.not.(associated(nodes))) then; get_next_=next_tag(cur_tag); return; end if
    nodes_index=nodes_index+1
    if (nodes_index>size(nodes)) then; get_next_=-1; return; end if
    cur_tag=>nodes(nodes_index)%node
    get_next_=size(nodes)-nodes_index
    
  end function

  function next_tag(cur_tag)
    type(tag),pointer :: cur_tag
    integer,save :: level,prev_index
    integer :: next_tag
    
    if(associated(cur_tag%children)) then
       if(size(cur_tag%children)>prev_index) then;          
          cur_tag=>cur_tag%children(prev_index+1)
          prev_index=1
          next_tag=1
          level=level+1
          return
       end if
    end if
    do while(associated(cur_tag%parent).and.level>0)
       prev_index=cur_tag%child
       cur_tag=>cur_tag%parent
       level=level-1
       if(prev_index<size(cur_tag%children)) then
          cur_tag=>cur_tag%children(prev_index+1)
          prev_index=1
          next_tag=1
          level=level+1
          return
       end if
    end do
    cur_tag=>null()
  end function


  integer function get_next_named(cur_tag,name,id,value)
    type(tag),pointer :: cur_tag
    character(len=*) :: name
    character(len=*),optional :: id,value
    integer :: i
    integer,save :: cur_level,prev_index
    type(tag_ref),pointer :: nodes(:)
    
    if(name=='') then;
       get_next_named=next_tag(cur_tag)
       return
    end if
    nodes=>get_node(cur_tag,name,id,value)
    nodes_index=0
    if(size(nodes)>0) cur_tag=>nodes(nodes_index)%node
    get_next_named=size(nodes)-nodes_index
  end function



  function xml_close(doc_name)
     character(len=*),optional :: doc_name
     type(tag),pointer         :: xml_close
     type(tag),pointer         :: tmp_document
     type(tag),pointer         :: tmp_xml_documents(:)
     integer :: i
     
     if (present(doc_name)) then
        do i=1,size(xml_documents)
           if (xml_documents(i)%name==doc_name) exit
        end do
        if (.not.(xml_documents(i)%name==doc_name)) then
           write(*,*) 'Warning: document name not a valid XML tree: ',doc_name
           return
        end if
     else; i=size(xml_documents) 
     end if
     tmp_document=>xml_documents(i)
     write(*,*) 'Deleting children'
     tmp_document=>del_children(tmp_document)
     write(*,*) 'Deleting documents'
     if (size(xml_documents)>1) then
        write(*,*) 'tst2'

	if (present(doc_name)) then
           do i=1,size(xml_documents)	
              if(doc_name==xml_documents(i)%name) exit
           end do
        else; i=size(xml_documents)
        end if
        allocate(tmp_xml_documents(size(xml_documents)-1))
        tmp_xml_documents(1:i-1)=xml_documents(1:i-1)
        tmp_xml_documents(i:)=xml_documents(i+1:)
        xml_documents=>tmp_xml_documents
        write(*,*) 'tst2'
     else
        write(*,*) 'tst'
!        deallocate(xml_documents(1))
        xml_documents=>null()
        write(*,*) 'tst'
     end if
     xml_close=>null()
  end function

 
  function xml_open(doc_name)
    type(tag),pointer :: xml_open
    character(len=*)  :: doc_name
    type(tag),pointer :: tmp_document(:)
    integer           :: new_size

! Instantiate a new XML document tree.
    new_size=1
    if (associated(xml_documents)) then
       new_size=size(xml_documents)+1       
       allocate(tmp_document(new_size))
       tmp_document=xml_documents
       deallocate(xml_documents)
       xml_documents=>tmp_document
    else; allocate(xml_documents(new_size))
    end if        
! Name the document
    
    xml_documents(new_size)%name=doc_name
    xml_open=>null()
    xml_open=>xml_documents(new_size)
    nullify(xml_open%data,xml_open%real_data,xml_open%attr,xml_open%children)
  end function

  function xml_read(file_name,doc_name)
    use xml_file

    character(len=*)         :: file_name
    character(len=*),optional:: doc_name
    character(len=name_len)  :: var
    character(len=attr_len)  :: attr
    character(len=1)         :: chr
    type(tag),pointer        :: xml_read
    type(tag),pointer        :: node
    integer                  :: ios=0,file_num,i,j,align
    integer(kind=8)          :: data_size, data_width,new_size
    logical                  :: is_tag=.false.,new_node=.false.
    type(tag),pointer        :: tmp_xml_documents(:)
    
    if(present(doc_name)) then; xml_read=xml_open(doc_name)
    else; xml_read=xml_open(file_name)
    end if

    node=>xml_read
! Point to the newly created document tree

    file_num=open_file(file_name)

    if (file_num<0) return
    do while(ios==0.or.ios==-2)
       var=''
       ios=read_var(file_num,var)
       select case(var)
       case ('<')
          is_tag=.true.
          ios=read_var(file_num,var)
          new_node=.true.
       case ('>')
          is_tag=.false.
          if (node%name=='dataset' .or. get_attr(node,'type')=='data') then
          ! Read a complet chunk of real numbers.
	     attr=get_attr(node,'size')
             read(attr, '(i8)') data_size
	     attr=get_attr(node,'width')
             read(attr, '(i8)') data_width
             if(data_width==0) data_width=1
             allocate(node%real_data(data_width,data_size))
             ! File pointer to next line
             read(file_num,'(a)',advance='yes') var

!            determine tab level
             chr=' '
             i=-1
             do while (chr==' ')
                read(file_num,'(1a)',advance='no') chr
                i=i+1
             end do
             backspace(file_num)

             if (i>0.and.(.not.(chr=='-'))) then; align=i-1
             else; align=i; end if             

             do i=1,data_size
                do j=1,align
                   read(file_num,'(a)',advance='no') chr
                end do
                do j=1,data_width-1
                   read(file_num,real_fmt,advance='no') node%real_data(j,i)
                end do
                read(file_num,real_fmt,advance='yes') node%real_data(data_width,i)
             end do
          end if
        case(' ')
           i=0
        case default
          if (is_tag) then          
             ios=read_var(file_num,attr)
             if (.not.(attr)=='=') then 
                write(*,*) 'ERROR: attribute not of the form: name="value"'
                write(*,*) 'TAG:      ', node%name
                write(*,*) 'ATTRIBUTE:', var
             end if
             ios=read_var(file_num,attr)
             call add_attr(node,var,attr)        
          else
             if (.not.(trim(var))=='') call add_data(node,trim(var))
          end if
       end select

       if (new_node) then
          if (var(1:1)=='/') then
            if (var(2:)==node%name(:)) then
               if (associated(node%parent)) then
                  node=>node%parent
                end if
            else
               write(*,*) "ERROR: Opening and closing tags are different"
               write(*,*) "Opening tag: ", node%name
               write(*,*) "Closing tag: ", var
               node=>xml_documents(size(xml_documents))
               node=>del_children(node)               
               return
            end if
          else
             node=>add_child(node,var)         
          end if
          new_node=.false.
       end if
    end do
    close(file_num)
    write(*,*) 'Finished reading'
  end function
  
  recursive subroutine write_child(file_num,node)
    use xml_file

    type(tag),pointer :: node,tmp_node
    integer :: file_num
    integer,save :: tab_level=0
    character(len=attr_len) :: fmt

    integer :: i,j



    if(tab_level>0) then;
      write(fmt,'("(t",i4,",a)")') tab_level
    else; fmt='(a)'; end if;

    write(file_num,fmt,advance='no') '<'//trim(node%name)
    if(associated(node%attr)) then
       do i=1,size(node%attr,2)
          write(file_num, '(a)',advance='no') ' '//trim(node%attr(1,i))//'="'//trim(node%attr(2,i))//'"'
       end do
    end if
    write(file_num,'(a)',advance='yes') '>'
    if (associated(node%children)) then
       tab_level=tab_level+tab_default
       do i=1,size(node%children)
          tmp_node=>node%children(i)
          call write_child(file_num,tmp_node)
       end do
       tab_level=tab_level+tab_default
    end if
    if (associated(node%real_data)) then
       do i=1,size(node%real_data,2) 
          write(file_num,fmt,advance='no')
          do j=1,size(node%real_data,1)-1
             write(file_num,real_fmt,advance='no') node%real_data(j,i)
          end do
          write(file_num,real_fmt,advance='yes') node%real_data(j,i)
       end do
    end if
    if (associated(node%data)) then
       write(file_num,fmt,advance='no')
       do i=1,size(node%data)
          write(file_num,'(a)',advance='no') node%data(i)
       end do
    end if
    write(file_num,fmt,advance='yes') '</'//trim(node%name)//'>'
    

  end subroutine


  subroutine xml_write(file_name,doc_name)
    use xml_file

    character(len=*)  :: file_name
    character(len=*),optional  :: doc_name
    integer :: file_num, doc_num,i
    type(tag), pointer :: tmp_doc,tmp_node
    type(tag),pointer :: tmp
    
  
    file_num=open_file(file_name)
    if (file_num<0) then;
       write(*,*) 'Failed to open file.'
       return
    end if
    if (present(doc_name)) then
       do doc_num=1,size(xml_documents)
          write(*,*) 'name:', xml_documents(doc_num)%name
          if (xml_documents(doc_num)%name==doc_name) exit
       end do
       if (.not.(xml_documents(doc_num)%name==doc_name)) then
          write(*,*) 'No document with name: ', doc_name
          return
       end if
    else; doc_num=size(xml_documents)
    end if 
    write(*,*) 'doc_num', doc_num
    tmp_doc=>xml_documents(doc_num)
    do i=1,size(tmp_doc%children)
       tmp_node=>tmp_doc%children(i)
       call write_child(file_num,tmp_node)
    end do
    close(file_num)

  end subroutine
  


end module
