module Taiji.Core.Utils
    ( SiteAffinity
    , toSiteAffinity
    , getSiteAffinity
    , PeakAffinity
    , toPeakAffinity
    , getPeakAffinity
    , SiteInfo(..)
    , TFBS
    , readExpression
    , lp
    ) where

import qualified Data.HashMap.Strict     as M
import qualified Data.ByteString.Char8   as B
import qualified Data.Vector.Unboxed as U
import           Data.CaseInsensitive              (CI)
import Bio.Data.Bed
import Data.List (foldl')
import           Data.CaseInsensitive    (mk)
import           Data.Maybe              (fromMaybe)
import           Bio.Utils.Misc          (readDouble)
import           Bio.Utils.Functions               (scale)
import qualified Data.Matrix.Unboxed               as MU

import           Taiji.Types

-- | Affinity score of a TF binding site, from 0 to 1.
newtype SiteAffinity = SiteAffinity
    { getSiteAffinity :: Double } deriving (Ord, Eq)

-- | Convert score [0,1000] to affinity score [0,1].
toSiteAffinity :: Int -> SiteAffinity
toSiteAffinity x = SiteAffinity $ fromIntegral x / 1000
{-# INLINE toSiteAffinity #-}

-- | Affinity score of a peak, from 0 to 1.
newtype PeakAffinity = PeakAffinity
    { getPeakAffinity :: Double } deriving (Ord, Eq)

-- | Convert p-value to affinity score [0,1].
toPeakAffinity :: Double -> PeakAffinity
toPeakAffinity x = PeakAffinity $ 1 / (1 + exp (-(x - 5)))
{-# INLINE toPeakAffinity #-}

data SiteInfo = SiteInfo
    { _tf_name :: CI B.ByteString
    , _site_affinity :: SiteAffinity
    , _peak_affinity :: PeakAffinity }

type TFBS = BEDExt BED3 SiteInfo


-- | Read RNA expression data
readExpression :: Double    -- ^ Threshold to call a gene as non-expressed
               -> B.ByteString  -- ^ cell type
               -> FilePath
               -> IO (M.HashMap GeneName (Double, Double)) -- ^ absolute value and z-score
readExpression cutoff ct fl = do
    c <- B.readFile fl
    let ((_:header):dat) = map (B.split '\t') $ B.lines c
        rowNames = map (mk . head) dat
        dataTable = map (U.fromList . map ((+pseudoCount) . readDouble) . tail) dat
        dataTable' = MU.fromRows $ zipWith U.zip dataTable $
            map computeZscore dataTable :: MU.Matrix (Double, Double)
        idx = fromMaybe (error $ "Cell type:" ++ B.unpack ct ++ " not found!") $
            lookup ct $ zip header [0..]
    return $ M.fromList $ zip rowNames $ U.toList $ dataTable' `MU.takeColumn` idx
  where
    pseudoCount = 0.1
    computeZscore xs
        | U.length xs == 1 = U.map log xs
        | U.all (<cutoff) xs = U.replicate (U.length xs) 0
        | U.length xs == 2 = let fc = log $ U.head xs / U.last xs
                             in U.fromList [fc, negate fc]
        | U.all (== U.head xs) xs = U.replicate (U.length xs) 0
        | otherwise = scale xs
{-# INLINE readExpression #-}

lp :: Int -> [Double] -> Double
lp p = (**(1/fromIntegral p)) . foldl' (+) 0 . map (**fromIntegral p)
{-# INLINE lp #-}