namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// Window内の変異数集計情報
    /// </summary>
    internal class VariantCount
    {
        /// <summary>
        /// Window内の変異数集計情報
        /// </summary>
        /// <param name="count">window内の変異数</param>
        /// <param name="bulk1ZeroSnpIndexVariantCount">Bulk1 SNP-indexが0の変異数</param>
        /// <param name="bulk2ZeroSnpIndexVariantCount">Bulk2 SNP-indexが0の変異数</param>
        public VariantCount(int count, int bulk1ZeroSnpIndexVariantCount, int bulk2ZeroSnpIndexVariantCount)
        {
            Count = count;
            Bulk1ZeroSnpIndexVariantCount = bulk1ZeroSnpIndexVariantCount;
            Bulk2ZeroSnpIndexVariantCount = bulk2ZeroSnpIndexVariantCount;
        }


        /// <summary>
        /// 変異数を取得する。
        /// </summary>
        public int Count { get; }

        /// <summary>
        /// Bulk1でSNP-indexが0になる変異数を取得する。
        /// </summary>
        public int Bulk1ZeroSnpIndexVariantCount { get; }

        /// <summary>
        /// Bulk2でSNP-indexが0になる変異数を取得する。
        /// </summary>
        public int Bulk2ZeroSnpIndexVariantCount { get; }
    }
}
