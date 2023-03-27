namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// RefAltタイプリードカウンター
    /// </summary>
    internal class RefAltReadCounter : ReadCounter
    {
        /// <summary>
        /// リード数をカウントする。
        /// </summary>
        /// <param name="ad">AD</param>
        /// <returns>(refCount, altCount)</returns>
        public override (int refCount, int altCount) Count(AD ad)
        {
            return (ad.ReadCounts[0], ad.ReadCounts[1]);
        }
    }
}
