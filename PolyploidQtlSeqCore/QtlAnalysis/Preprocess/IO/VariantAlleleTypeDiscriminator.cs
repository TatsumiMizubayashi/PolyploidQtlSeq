namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// 変異アレルタイプ識別機
    /// </summary>
    internal static class VariantAlleleTypeDiscriminator
    {
        /// <summary>
        /// 変異アレルタイプを識別する。
        /// </summary>
        /// <param name="altAlleles">Altアレル配列</param>
        /// <param name="ads">AD配列</param>
        /// <returns>アレルの種類</returns>
        public static VariantAlleleType Discern(string[] altAlleles, params AD[] ads)
        {
            if (altAlleles.Length == 1) return VariantAlleleType.RefAlt;
            if (altAlleles.Length >= 3) return VariantAlleleType.Multi;

            // アレルが3個の場合
            var countDictionary = CreateCountDictionary(ads);
            var hasReadAlleleCount = countDictionary.Values.Count(x => x > 0);
            return hasReadAlleleCount == 2  // REFが0個の場合
                ? VariantAlleleType.Alt1Alt2
                : VariantAlleleType.Multi;
        }

        private static IReadOnlyDictionary<int, int> CreateCountDictionary(AD[] ads)
        {
            // アレルは３個限定
            var countDictionary = new Dictionary<int, int>()
            {
                [0] = 0,
                [1] = 0,
                [2] = 0
            };

            foreach (var ad in ads)
            {
                if (ad.IsNoData) continue;

                for (var i = 0; i < ad.ReadCounts.Length; i++)
                {
                    countDictionary[i] += ad.ReadCounts[i];
                }
            }

            return countDictionary;
        }
    }
}
