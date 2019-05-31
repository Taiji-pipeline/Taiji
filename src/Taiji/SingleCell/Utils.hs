module Taiji.SingleCell.Utils where

import Bio.Data.Bed
import qualified Data.HashMap.Strict                   as M

import Taiji.Prelude

mkTFBSMap :: Monad m => ConduitT TFBS o m (BEDTree [SiteInfo])
mkTFBSMap = mapC (\x -> (x^._bed, [x^._data])) .| sinkList >>=
    return . (fmap . fmap) nub' . bedToTree (++)
  where
    nub' = M.elems . M.fromListWith (\a b -> if g a > g b then a else b) .
        map (\x -> (_tf_name x, x))
      where
        g = getSiteAffinity . _site_affinity
{-# INLINE mkTFBSMap #-}