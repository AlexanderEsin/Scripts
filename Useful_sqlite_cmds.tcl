select DISTINCT binomial from t1 where binomial GLOB '*__*';
UPDATE t1 SET binomial = $sdash_name WHERE binomial = $ddash_name
select * from t1 where binomial = 'Alistipes_shahii' LIMIT 1;