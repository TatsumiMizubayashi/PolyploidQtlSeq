namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// AltAltタイプのリードカウンター
    /// </summary>
    internal class AltAltReadCounter : ReadCounter
    {
        /*
         * このタイプの場合はALTアレルが2種類あり、Ref型アレルが検出されない
         * 
         */

        /// <summary>
        /// リード数をカウントする。
        /// </summary>
        /// <param name="ad">AD</param>
        /// <returns>(alt1Count, alt2Count)</returns>
        public override (int refCount, int altCount) Count(AD ad)
        {
            return (ad.ReadCounts[1], ad.ReadCounts[2]);
        }
    }
}
