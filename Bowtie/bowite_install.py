# 기존 파일들 삭제
rm -f GRCh38_noalt_as.zip*

# 새로 다운로드
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

# 파일 확인
ls -lh GRCh38_noalt_as.zip

# 압축 해제
unzip GRCh38_noalt_as.zip
